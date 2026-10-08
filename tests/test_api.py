"""HTTP-level tests: both shapes, error taxonomy, conflicts, caching."""
import logging

from fastapi.testclient import TestClient

import main

logger = logging.getLogger(__name__)
client = TestClient(main.app)


def test_health_and_root():
    assert client.get("/health").json()["offline"] is True
    assert client.get("/").json()["service"]


def test_normalize_single_sample_ok():
    resp = client.post("/normalize", json={"text": "chr1 100 A T 0/1", "sample": "S1"})
    assert resp.status_code == 200
    body = resp.json()
    logger.info("输入摘要 /normalize 单行 -> valid=%s entry=%s",
                body["valid"], body["entries"][0]["variant"])
    assert body["valid"] is True
    assert body["entries"][0]["samples"] == ["S1"]
    assert body["entries"][0]["genotype_by_sample"]["S1"] == "0/1"


def test_normalize_malformed_reports_distinguishable_reason_and_location():
    text = "chr1 x A T 0/1\nchr1 100 A U 0/1\nchr1 100 A T 9/9\nshort"
    resp = client.post("/normalize", json={"text": text, "sample": "S1"})
    body = resp.json()
    reasons = [(r["reason"], r["location"]["line_no"]) for r in body["rejected"]]
    logger.info("输入摘要 4 行畸形 -> %s", reasons)
    assert body["valid"] is False
    assert [r[0] for r in reasons] == ["BAD_POSITION", "BAD_ALLELE", "BAD_GENOTYPE", "MISSING_FIELDS"]
    assert [r[1] for r in reasons] == [1, 2, 3, 4]
    assert all(r["location"]["block"] == "sample:S1" for r in body["rejected"])


def test_normalize_equivalent_writings_converge():
    text = "1 100 ACG AT 0/1\nchr1 101 CG T 1/1"
    body = client.post("/normalize", json={"text": text, "sample": "S1"}).json()
    logger.info("等价写法 -> entries=%d", len(body["entries"]))
    assert len(body["entries"]) == 1
    assert body["entries"][0]["variant"] == "CHR1:101:CG:T"
    assert body["entries"][0]["conflicts"][0]["type"] == "GENOTYPE_CONFLICT"


def test_callset_multi_sample_aggregation_and_conflict():
    text = (
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\t0/1\n"          # two carriers, one entry
        "chr1\t100\tA\tT\t1/1\t0/1\n"          # S1 HET vs HOM -> conflict
        "bad row here\n"
    )
    body = client.post("/callset", json={"text": text}).json()
    logger.info("输入摘要 多样本 -> entries=%d conflicts=%s rejected=%d",
                len(body["entries"]), [c["type"] for c in body["conflicts"]],
                len(body["rejected"]))
    assert len(body["entries"]) == 1
    assert body["entries"][0]["samples"] == ["S1", "S2"]
    assert [c["type"] for c in body["conflicts"]] == ["GENOTYPE_CONFLICT"]
    assert body["conflicts"][0]["sample"] == "S1"
    assert body["rejected"][0]["reason"] == "MISSING_FIELDS"


def test_callset_ref_conflict_marked():
    text = (
        "#CHROM\tPOS\tREF\tALT\tS1\n"
        "chr1\t100\tA\tT\t1/1\n"
        "chr1\t100\tA\tT\t0/0\n"
    )
    body = client.post("/callset", json={"text": text}).json()
    assert {c["type"] for c in body["conflicts"]} == {"REF_CONFLICT"}
    assert body["valid"] is False


def test_empty_payload_is_422():
    assert client.post("/normalize", json={"text": "  ", "sample": "S1"}).status_code == 422
    assert client.post("/normalize", json={"text": "x", "sample": " "}).status_code == 422
    assert client.post("/callset", json={"text": ""}).status_code == 422


def test_repeated_post_is_stable_and_cached():
    text = "chr1 100 A T 0/1\nchr2 200 G C 1/1"
    bodies = [client.post("/normalize", json={"text": text, "sample": "S9"}).json()
              for _ in range(3)]
    sig = [(e["variant"], e["samples"]) for e in bodies[0]["entries"]]
    for b in bodies:
        assert [(e["variant"], e["samples"]) for e in b["entries"]] == sig
    assert bodies[0]["content_hash"] == bodies[-1]["content_hash"]


def test_callset_accepts_single_sample_shaped_text():
    # single-sample text (4 cols) via /callset -> one inferred sample
    body = client.post("/callset", json={"text": "chr1 100 A T 0/1"}).json()
    assert body["entries"][0]["genotype_by_sample"] == {"sample_1": "0/1"}


def _baseline():
    body = client.post("/callset", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS1\tS2\nchr1\t100\tA\tT\t0/1\t1/1\n",
    }).json()
    return body["version_id"], body


def test_versioned_callset_exposes_version_and_initial_summary():
    _, body = _baseline()
    logger.info("判定依据 初始版本: version_id 存在=%s, parent=%s, summary.baseline=%s",
                bool(body["version_id"]), body["parent_version_id"],
                body["change_summary"]["baseline_version"])
    assert len(body["version_id"]) == 64
    assert body["parent_version_id"] is None
    assert body["change_summary"]["baseline_version"] is None
    assert body["change_summary"]["entries_added"] == ["CHR1:100:A:T"]


def test_increment_endpoint_multi_and_single_shape_diff():
    vid, _ = _baseline()
    body = client.post("/increment", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS3\tS1\nchr1\t100\tA\tT\t1/1\t1/1\n"
                "chr3\t300\tC\tG\t0/1\t0/1\n",
        "baseline_version": vid,
    }).json()
    cs = body["change_summary"]
    logger.info("判定依据 /increment 多样本: added=%s modified变体=%s new=%s corrected=%s",
                cs["entries_added"], [m["variant"] for m in cs["entries_modified"]],
                cs["new_samples"], cs["corrected_samples"])
    assert body["parent_version_id"] == vid
    assert cs["entries_added"] == ["CHR3:300:C:G"]
    assert cs["new_samples"] == ["S3"] and cs["corrected_samples"] == ["S1"]
    modified = {m["variant"]: m for m in cs["entries_modified"]}
    changes = {(x["change"], x["sample"]): (x["before"], x["after"])
               for x in modified["CHR1:100:A:T"]["sample_changes"]}
    assert changes[("GENOTYPE_CHANGED", "S1")] == ("0/1", "1/1")
    assert changes[("CARRIER_ADDED", "S3")] == (None, "1/1")

    # single-sample shape correction through the same endpoint
    body2 = client.post("/increment", json={
        "text": "chr1 100 A T 0/1", "sample": "S1",
        "baseline_version": body["version_id"],
    }).json()
    chr1 = body2["entries"][0]
    assert chr1["genotype_by_sample"]["S1"] == "0/1"
    logger.info("判定依据 /increment 单样本形态修正 S1 -> 0/1，未与旧 1/1 叠加")


def test_increment_idempotent_over_http():
    vid, _ = _baseline()
    payload = {
        "text": "#CHROM\tPOS\tREF\tALT\tS3\nchr3\t300\tC\tG\t0/1\n",
        "baseline_version": vid,
    }
    responses = [client.post("/increment", json=payload).json() for _ in range(3)]
    versions = {b["version_id"] for b in responses}
    summaries = [b["change_summary"] for b in responses]
    logger.info("判定依据 重复增量: 唯一版本数=%d, 摘要一致=%s",
                len(versions), all(s == summaries[0] for s in summaries))
    assert versions == {responses[0]["version_id"]}
    assert all(s == summaries[0] for s in summaries)


def test_increment_error_codes_are_distinguishable():
    vid, _ = _baseline()
    r404 = client.post("/increment", json={
        "text": "chr1 100 A T 0/1", "sample": "S1", "baseline_version": "missing"})
    assert r404.status_code == 404
    assert r404.json()["detail"]["code"] == "BASELINE_NOT_FOUND"

    r422 = client.post("/increment", json={"text": "# c\n", "baseline_version": vid})
    assert r422.status_code == 422
    assert r422.json()["detail"]["code"] == "EMPTY_INCREMENT"
    logger.info("判定依据 /increment 错误分类: 404 BASELINE_NOT_FOUND / 422 EMPTY_INCREMENT")

    # malformed delta stays a normal 200 valid=false response, not swallowed
    bad = client.post("/increment", json={
        "text": "#CHROM\tPOS\tREF\tALT\tS9\nchr9\t900\tA\tT\t0/1\n"
                "chr8\t800\tA\tT\t0/1\tEXTRA\n",
        "baseline_version": vid,
    }).json()
    assert bad["valid"] is False
    assert [r["reason"] for r in bad["rejected"]] == ["HEADER_SAMPLE_MISMATCH"]
    logger.info("判定依据 /increment 列错位 -> %s（不被静默吃入）",
                [r["reason"] for r in bad["rejected"]])
