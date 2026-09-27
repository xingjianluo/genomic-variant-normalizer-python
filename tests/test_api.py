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
