"""api 端到端测试：双形态输入、拒绝分类、冲突标记、重复提交一致性、并发。"""
import concurrent.futures
import logging

from fastapi.testclient import TestClient

from main import app

logger = logging.getLogger("test.api")
client = TestClient(app)


def test_health_is_offline():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["offline"] is True


def test_single_and_vcf_forms_unify():
    single = [{
        "sample_id": "S1",
        "format": "single",
        "text": "chr1\t10\tA\tG\t0/1",
    }]
    vcf = [{
        "sample_id": "",
        "format": "vcf",
        "text": "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1",
    }]
    r1 = client.post("/callsets/inspect", json=single).json()
    r2 = client.post("/callsets/inspect", json=vcf).json()
    k1 = (r1["variants"][0]["chrom"], r1["variants"][0]["pos"],
          r1["variants"][0]["ref"], r1["variants"][0]["alt"])
    k2 = (r2["variants"][0]["chrom"], r2["variants"][0]["pos"],
          r2["variants"][0]["ref"], r2["variants"][0]["alt"])
    logger.info("输入摘要 single/vcf 两种形态 -> 归一键 %s / %s，一致=%s", k1, k2, k1 == k2)
    assert k1 == k2


def test_malformed_input_returns_distinguishable_rejections():
    payload = [{
        "sample_id": "S1",
        "format": "single",
        "text": "\n".join([
            "chr1\t100\tA",            # 列数错
            "chr1\tx\tA\tG",           # 位置错
            "chrZ\t1\tA\tG",           # 染色体错
            "chr1\t1\tA\t<DEL>",       # 等位错
        ]),
    }]
    resp = client.post("/callsets/inspect", json=payload)
    assert resp.status_code == 200  # 畸形记录不导致整体失败
    data = resp.json()
    reasons = {r["reason"] for r in data["rejected_records"]}
    logger.info("输入摘要 4 条畸形行 -> 可区分原因=%s", reasons)
    assert reasons == {"WRONG_COLUMN_COUNT", "INVALID_POSITION",
                       "INVALID_CHROM", "INVALID_ALLELE"}
    assert all(r["line"] >= 1 and r["detail"] for r in data["rejected_records"])
    summary = client.post("/callsets", json=payload).json()
    assert summary["rejected"] == 4 and summary["accepted_calls"] == 0


def test_conflict_is_marked_end_to_end():
    payload = [{"sample_id": "S1", "format": "single", "text": "\n".join([
        "chr1\t10\tA\tG\tHET",
        "chr1\t10\tA\tG\t1/1",
    ])}]
    data = client.post("/callsets/inspect", json=payload).json()
    v = data["variants"][0]
    logger.info("输入摘要 同样本 HET+HOM_ALT -> conflicting=%s conflicts=%s",
                v["conflicting"], v["conflicts"])
    assert v["conflicting"] is True
    assert v["conflict_types"] == ["GENOTYPE_CONFLICT"]
    assert set(v["genotype_summary"]["S1"]) == {"HET", "HOM_ALT"}


def test_duplicate_submission_is_stable_and_cached():
    payload = [{"sample_id": "S1", "format": "single",
                "text": "chr1\t10\tA\tG\tHET\nchr1\t10\tA\tG\tHET"}]
    a = client.post("/callsets/inspect", json=payload).json()
    b = client.post("/callsets/inspect", json=payload).json()
    assert a == b
    assert a["accepted_calls"] == 1
    fetched = client.get(f"/callsets/{a['callset_id']}").json()
    assert fetched == a
    logger.info("判定: 重复提交响应完全一致, callset_id=%s 可 GET 复用", a["callset_id"])
    assert client.get("/callsets/nonexistent").status_code == 404


def test_empty_body_rejected():
    assert client.post("/callsets", json=[]).status_code == 400


def test_concurrent_requests_no_pollution():
    def post(sample, pos):
        body = [{"sample_id": sample, "format": "single",
                 "text": f"chr1\t{pos}\tA\tG\tHET"}]
        return client.post("/callsets/inspect", json=body).json()

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(lambda i: post(f"S{i % 3}", 100 + i), range(24)))
    ids = [r["callset_id"] for r in results]
    # 相同 (sample,pos) 的请求结论一致；不同输入结论不同
    groups: dict[str, set] = {}
    for i, cid in enumerate(ids):
        groups.setdefault((f"S{i % 3}", 100 + i), set()).add(cid)
    assert all(len(v) == 1 for v in groups.values())
    for r in results:
        samples = {s for v in r["variants"] for s in v["samples"]}
        assert len(samples) == 1
    logger.info("输入摘要 24 个并发请求(3 样本轮转) -> 唯一 callset 数=%d，无污染",
                len(set(ids)))
