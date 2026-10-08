"""HTTP tests for incremental endpoints: both shapes, versions, errors.

Endpoints:
* ``POST /increment/root/single``  {text, sample}            -> version 0
* ``POST /increment/root/multi``   {text}                    -> version 0
* ``POST /increment/single``       {text, sample, base_version}
* ``POST /increment/multi``        {text, base_version}
* ``GET  /version/{id}``
"""
import logging
from concurrent.futures import ThreadPoolExecutor

from fastapi.testclient import TestClient

import main
from variant_service.api import service

logger = logging.getLogger(__name__)
client = TestClient(main.app)


def test_root_then_single_increment_preserves_untouched():
    root = client.post("/increment/root/multi", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS1\tS2\nchr1\t100\tA\tT\t0/1\t1/1\n"
    }).json()
    resp = client.post("/increment/single", json={
        "text": "chr3 300 C G 0/1", "sample": "S3",
        "base_version": root["version_id"],
    })
    assert resp.status_code == 200
    body = resp.json()
    logger.info(
        "判定依据[/increment/single 新增样本]: version=%s 前 8 位 seq=%d added=%s samples_added=%s",
        body["version_id"][:8], body["sequence"],
        body["changes"]["added_entries"], body["changes"]["samples_added"],
    )
    assert body["version_id"] != root["version_id"]
    assert body["parent_version"] == root["version_id"]
    assert body["sequence"] == 1
    assert body["changes"]["added_entries"] == ["CHR3:300:C:G"]
    assert body["changes"]["samples_added"] == ["S3"]
    first = body["entries"][0]
    assert first["variant"] == "CHR1:100:A:T"
    assert first["samples"] == ["S1", "S2"]


def test_increment_multi_shape_and_version_get():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "MS1"}).json()
    text = (
        "#CHROM\tPOS\tREF\tALT\tMS2\tMS3\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
        "chr4\t400\tT\tA\t0/1\t./.\n"
    )
    body = client.post("/increment/multi",
                       json={"text": text, "base_version": root["version_id"]}).json()
    logger.info("判定依据[/increment/multi 多样本]: entries=%d changed=%s added=%s",
                len(body["entries"]),
                [c["variant"] for c in body["changes"]["changed_entries"]],
                body["changes"]["added_entries"])
    changed = {c["variant"]: c for c in body["changes"]["changed_entries"]}
    assert changed["CHR1:100:A:T"]["samples_added"] == ["MS2", "MS3"]
    assert changed["CHR1:100:A:T"]["kind"] == "CHANGED"

    got = client.get(f"/version/{body['version_id']}")
    assert got.status_code == 200
    assert got.json()["version_id"] == body["version_id"]
    assert got.json()["entries"] == body["entries"]
    assert client.get("/version/unknown-id").status_code == 404


def test_correction_vs_evidence_stacking_over_http():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "CS1"}).json()
    body = client.post("/increment/single", json={
        "text": "chr1 100 A T 1/1", "sample": "CS1",
        "base_version": root["version_id"],
    }).json()
    gt = body["entries"][0]["genotype_by_sample"]["CS1"]
    conflict_types = [c["type"] for c in body["entries"][0]["conflicts"]]
    logger.info("判定依据[修正不叠加 HTTP]: CS1 修正后=%s, 条目内冲突=%s", gt, conflict_types)
    assert gt == "1/1"
    assert "GENOTYPE_CONFLICT" not in conflict_types
    change = body["changes"]["changed_entries"][0]
    assert change["genotype_changes"] == [{"sample": "CS1", "old": "0/1", "new": "1/1"}]


def test_malformed_increment_returns_distinguishable_reasons():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "BS1"}).json()
    body = client.post("/increment/single", json={
        "text": "chr1 x A T 0/1\nchr1 100 A N 0/1\nbad",
        "sample": "BS2", "base_version": root["version_id"],
    }).json()
    reasons = [(r["reason"], r["location"]["line_no"]) for r in body["rejected"]]
    logger.info("判定依据[畸形增量 HTTP 可区分]: %s, valid=%s", reasons, body["valid"])
    assert reasons == [("BAD_POSITION", 1), ("BAD_ALLELE", 2), ("MISSING_FIELDS", 3)]
    assert body["valid"] is False
    # old entry untouched
    assert body["entries"][0]["samples"] == ["BS1"]


def test_stale_baseline_is_409_with_reason_and_current_version():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "OL1"}).json()
    client.post("/increment/single", json={
        "text": "chr1 100 A T 1/1", "sample": "OL1",
        "base_version": root["version_id"],
    })
    resp = client.post("/increment/single", json={
        "text": "chr1 100 A G 0/1", "sample": "OL1",
        "base_version": root["version_id"],
    })
    assert resp.status_code == 409
    detail = resp.json()["detail"]
    logger.info("判定依据[过期基线 HTTP 409]: reason=%s current_version 存在=%s conflicting=%s",
                detail["reason"], bool(detail["current_version"]),
                detail["conflicting_samples"])
    assert detail["reason"] == "BASELINE_SUPERSEDED"
    assert detail["conflicting_samples"] == ["OL1"]
    assert detail["current_version"]


def test_unknown_baseline_is_404_and_empty_is_422():
    r1 = client.post("/increment/single", json={
        "text": "chr1 100 A T 0/1", "sample": "X1", "base_version": "nope",
    })
    assert r1.status_code == 404
    assert r1.json()["detail"]["reason"] == "BASELINE_NOT_FOUND"
    assert client.post("/increment/single", json={
        "text": "  ", "sample": "X1", "base_version": "x",
    }).status_code == 422
    assert client.post("/increment/multi", json={
        "text": " ", "base_version": "x",
    }).status_code == 422


def test_idempotent_replay_over_http():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "ID1"}).json()
    payload = {"text": "chr5 500 G C 0/1", "sample": "ID2",
               "base_version": root["version_id"]}
    a = client.post("/increment/single", json=payload).json()
    b = client.post("/increment/single", json=payload).json()
    logger.info("判定依据[HTTP 幂等]: 版本一致=%s, 摘要一致=%s",
                a["version_id"] == b["version_id"], a["changes"] == b["changes"])
    assert a["version_id"] == b["version_id"]
    assert a["changes"] == b["changes"]
    assert a["entries"] == b["entries"]


def test_empty_multi_increment_without_sample_columns_is_422():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "E1"}).json()
    # non-whitespace, non-comment line with the 4 fixed columns but no sample
    # column: nothing can be added or replaced -> EMPTY_INCREMENT
    resp = client.post("/increment/multi", json={
        "text": "chr1\t100\tA\tT\n",
        "base_version": root["version_id"],
    })
    assert resp.status_code == 422
    assert resp.json()["detail"]["reason"] == "EMPTY_INCREMENT"


def test_legacy_callset_baseline_can_be_extended():
    """The original endpoints return a version_id usable as increment baseline."""
    base = client.post("/callset", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS1\nchr1\t100\tA\tT\t0/1\n"
    }).json()
    assert base["version_id"]
    inc = client.post("/increment/multi", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS2\nchr1\t100\tA\tT\t1/1\n",
        "base_version": base["version_id"],
    }).json()
    logger.info("判定依据[老端点基线可增量]: %s -> %s, samples=%s",
                base["version_id"][:8], inc["version_id"][:8], inc["samples"])
    assert inc["parent_version"] == base["version_id"]
    assert inc["entries"][0]["samples"] == ["S1", "S2"]


def test_concurrent_overlapping_increments_http():
    root = client.post("/increment/root/single",
                       json={"text": "chr1 100 A T 0/1", "sample": "CC1"}).json()

    def job(i):
        r = client.post("/increment/single", json={
            "text": f"chr1 100 A T {'1/1' if i % 2 else '0/1'}",
            "sample": "CC1", "base_version": root["version_id"],
        })
        return r.status_code, (r.json().get("detail") or {}).get("reason")

    with ThreadPoolExecutor(max_workers=16) as pool:
        outcomes = list(pool.map(job, range(24)))
    rejects = [o for o in outcomes if o[0] == 409]
    oks = [o for o in outcomes if o[0] == 200]
    logger.info(
        "判定依据[HTTP 并发重叠]: 200=%d, 409=%d, 拒绝原因=%s (两种确定结论: "
        "命中幂等状态返回 200, 无法合并且状态未生成过时返回 409)",
        len(oks), len(rejects), {o[1] for o in rejects},
    )
    assert {o[0] for o in outcomes} <= {200, 409}
    assert {o[1] for o in rejects} <= {"BASELINE_SUPERSEDED"}
    assert oks  # at least one deterministic outcome always lands
    # head callset stays well-formed: exactly one sample, one entry
    root_snap = service.get_version(root["version_id"])
    head_id = service._lineage[root_snap.root_version]["head"]
    head = client.get(f"/version/{head_id}").json()
    assert head["samples"] == ["CC1"]
    assert len(head["entries"]) == 1
