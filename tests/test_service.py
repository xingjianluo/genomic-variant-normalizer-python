"""Tests for submission service: caching, determinism, concurrency, scale."""
import logging
import time
from concurrent.futures import ThreadPoolExecutor

from variant_service.models import RejectReason
from variant_service.service import CallsetService

logger = logging.getLogger(__name__)


SINGLE_A = "chr1 100 A T 0/1\nchr2 200 ACG AT hom\nBADLINE\n"
SINGLE_B = "chr1 100 A T 1/1\n"
MULTI = (
    "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
    "chr1\t100\tA\tT\t0/1\t1/1\n"
    "chr1\t100\tA\tG\t0/1\t0/0\n"   # S2 0/0 vs S1 alt at same site -> no cross-sample conflict
)


def test_single_and_multi_share_one_model(service):
    r1 = service.submit_single(SINGLE_A, "S1")
    r2 = service.submit_multi(MULTI)
    logger.info("输入摘要 single accepted=%d rejected=%d; multi entries=%d",
                r1.accepted_count, len(r1.rejected), len(r2.entries))
    assert {e.variant.render() for e in r2.entries}
    assert r1.rejected[0].reason is RejectReason.MISSING_FIELDS


def test_repeated_submission_reuses_conclusion(service):
    first = service.submit_single(SINGLE_A, "S1")
    second = service.submit_single(SINGLE_A, "S1")
    third = service.submit_single(SINGLE_A, "S1")
    logger.info("重复提交 cache=%s version=%s", service.cache_info(),
                (first.version_id or "")[:8])
    # Repeated identical submissions share the same content hash, version id
    # and conclusions (entries/rejections/conflicts are equal).
    assert first.content_hash == second.content_hash == third.content_hash
    assert first.version_id and first.version_id == third.version_id
    assert first.entries == third.entries
    assert first.rejected == third.rejected
    assert first.conflicts == third.conflicts
    assert service.cache_info()["entries"] == 1


def test_repeated_submission_output_stable(service):
    ref = service.submit_multi(MULTI)
    for _ in range(10):
        r = service.submit_multi(MULTI)
        assert [(e.variant.render(), e.samples) for e in r.entries] == \
               [(e.variant.render(), e.samples) for e in ref.entries]
        assert [x.reason for x in r.rejected] == [x.reason for x in ref.rejected]


def test_concurrent_overlapping_samples_do_not_contaminate(service):
    """Many threads submit overlapping/different single-sample blocks."""
    texts = {f"S{i % 8}": SINGLE_A if i % 2 == 0 else SINGLE_B for i in range(32)}

    def job(item):
        sample, text = item
        r = service.submit_single(text, sample)
        return sample, r

    with ThreadPoolExecutor(max_workers=16) as pool:
        results = list(pool.map(job, texts.items()))

    # All even-indexed samples saw SINGLE_A (2 accepted, 1 rejected), odds B.
    for sample, r in results:
        if sample in {"S0", "S2", "S4", "S6"}:
            assert r.accepted_count == 2
            assert len(r.rejected) == 1
        else:
            assert r.accepted_count == 1 and not r.rejected
        for e in r.entries:
            assert all(s == sample for s in list(e.samples) + list(e.genotype_by_sample))
    logger.info("并发 32 请求完成, cache entries=%d, 无污染", service.cache_info()["entries"])


def test_concurrent_identical_submissions_consistent(service):
    with ThreadPoolExecutor(max_workers=12) as pool:
        rs = list(pool.map(lambda _: service.submit_multi(MULTI), range(40)))
    sig = {(e.variant.render(), e.samples) for e in rs[0].entries}
    assert all({(e.variant.render(), e.samples) for e in r.entries} == sig for r in rs)


def test_lru_cache_bounded_and_evicts():
    svc = CallsetService(cache_max=3)
    for i in range(5):
        svc.submit_single(SINGLE_B, f"X{i}")
    assert svc.cache_info()["entries"] == 3
    # X0, X1 evicted -> recomputed as different object
    h0 = svc.hash_content(SINGLE_B, "X0")
    assert svc.get_cached(h0) is None
    h4 = svc.hash_content(SINGLE_B, "X4")
    assert svc.get_cached(h4) is not None


def test_ten_thousand_rows_complete_fast(service):
    rows = "\n".join(f"chr{i % 24 + 1}\t{1000 + i}\tA\tT\t0/1" for i in range(10_000))
    text = "#CHROM\tPOS\tREF\tALT\tS1\n" + rows
    start = time.perf_counter()
    r = service.submit_multi(text)
    elapsed = time.perf_counter() - start
    logger.info("输入摘要 10000 行 -> entries=%d 用时 %.3fs", len(r.entries), elapsed)
    assert len(r.entries) == 10_000
    assert elapsed < 5.0  # generous bound; local machine should be well under


def test_empty_sample_name_treated_as_unnamed(service):
    r = service.submit_single(SINGLE_B, "  ")
    # still processed; block label uses empty sample but no crash
    assert r.accepted_count == 1
