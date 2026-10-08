"""Service-level tests for incremental updates.

Business coverage (logged 判定依据):
* both input shapes as increments on a returned version_id;
* same baseline + same delta returns the same version object (idempotency);
* unknown baseline -> distinguishable BaselineNotFound; empty -> EmptyIncrement;
* concurrent overlapping increments: no sample loss, deterministic/reproducible;
* 10k-row scale: baseline build and increments stay fast and baseline reads are
  never blocked.
"""
import logging
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from variant_service.service import (
    BaselineNotFound,
    CallsetService,
    EmptyIncrement,
)

logger = logging.getLogger(__name__)

BASE = (
    "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
    "chr1\t100\tA\tT\t0/1\t1/1\n"
    "chr2\t200\tG\tC\t.\t0/0\n"
)
DELTA = (
    "#CHROM\tPOS\tREF\tALT\tS3\tS1\n"
    "chr1\t100\tA\tT\t1/1\t1/1\n"
    "chr3\t300\tC\tG\t0/1\t0/1\n"
)


@pytest.fixture
def svc():
    return CallsetService()


def test_initial_submissions_are_versioned(svc):
    v0 = svc.submit_multi_version(BASE)
    logger.info("判定依据 初始版本 version=%s parent=%s", v0.version_id[:10], v0.parent_version_id)
    assert v0.parent_version_id is None
    assert v0.summary.baseline_version is None
    assert set(v0.summary.entries_added) == {"CHR1:100:A:T", "CHR2:200:G:C"}
    # fetchable by id as a later baseline
    assert svc.get_version(v0.version_id) is v0


def test_multi_and_single_shape_increments(svc):
    v0 = svc.submit_multi_version(BASE)
    v1 = svc.increment_multi(DELTA, v0.version_id)
    logger.info("判定依据 多样本增量: added=%s modified=%d new=%s corrected=%s",
                v1.summary.entries_added, len(v1.summary.entry_changes),
                v1.summary.new_samples, v1.summary.corrected_samples)
    assert v1.parent_version_id == v0.version_id
    assert v1.summary.entries_added == ("CHR3:300:C:G",)
    assert v1.summary.new_samples == ("S3",)
    assert v1.summary.corrected_samples == ("S1",)

    v2 = svc.increment_single("chr1 100 A T 0/1", "S1", v1.version_id)
    chr1 = next(e for e in v2.result.entries if e.variant.pos == 100)
    logger.info("判定依据 单样本形态修正 S1 -> gt=%s（旧 HOM 不叠加）",
                chr1.genotype_by_sample["S1"].value)
    assert chr1.genotype_by_sample["S1"].value == "0/1"
    # S1's chr3 call from the multi delta is replaced away wholesale
    chr3 = next(e for e in v2.result.entries if e.variant.pos == 300)
    assert "S1" not in chr3.genotype_by_sample and "S3" in chr3.genotype_by_sample


def test_repeated_increment_is_idempotent_same_object(svc):
    v0 = svc.submit_multi_version(BASE)
    v1 = svc.increment_multi(DELTA, v0.version_id)
    v1b = svc.increment_multi(DELTA, v0.version_id)
    v1c = svc.increment_multi(DELTA, v0.version_id)
    logger.info("判定依据 同基线同增量重复提交: 同一对象=%s version=%s",
                v1 is v1b is v1c, v1.version_id[:10])
    assert v1 is v1b is v1c
    assert v1.summary == v1b.summary


def test_unknown_baseline_and_empty_increment_are_distinguishable(svc):
    v0 = svc.submit_multi_version(BASE)
    with pytest.raises(BaselineNotFound):
        svc.increment_multi("chr1 100 A T 0/1", "unknown-version-id")
    logger.info("判定依据 未知基线 -> BaselineNotFound(code=%s)", BaselineNotFound.code)
    with pytest.raises(EmptyIncrement):
        svc.increment_multi("# comment only\n", v0.version_id)
    with pytest.raises(EmptyIncrement):
        svc.increment_single("chr1 100 A T 0/1", "  ", v0.version_id)
    logger.info("判定依据 空增量/空样本 -> EmptyIncrement(code=%s)", EmptyIncrement.code)


def test_concurrent_overlapping_increments_lose_no_sample(svc):
    v0 = svc.submit_multi_version(BASE)

    def job(t):
        # half the threads correct an existing overlapping sample S1/S2;
        # half add brand-new samples; all branch from the SAME baseline.
        if t % 2 == 0:
            sample = ("S1", "S2")[t % 4 // 2]
            text = f"chr1 100 A T 1/1\nchr7 {700 + t} A T 0/1\n"
            return t, svc.increment_single(text, sample, v0.version_id)
        sample = f"N{t}"
        return t, svc.increment_single(f"chr8 {800 + t} A T 0/1\n", sample, v0.version_id)

    with ThreadPoolExecutor(max_workers=16) as pool:
        first = list(pool.map(job, range(40)))
    with ThreadPoolExecutor(max_workers=16) as pool:
        second = list(pool.map(job, range(40)))

    ids_first = [v.version_id for _, v in first]
    ids_second = [v.version_id for _, v in second]
    logger.info("判定依据 40 并发增量(样本重叠) 两次提交版本序列一致=%s",
                ids_first == ids_second)
    assert ids_first == ids_second  # reproducible under concurrency
    for t, v in first:
        present = {s for e in v.result.entries for s in e.genotype_by_sample}
        # baseline samples never disappear from a derived version
        assert {"S1", "S2"} <= present
        assert (f"N{t}" in present) if t % 2 else True
    # idempotent objects under concurrent duplicate submission
    same = sum(1 for (_, a), (_, b) in zip(first, second) if a is b)
    logger.info("判定依据 并发重复增量对象同一性复用 %d/40", same)
    assert same == 40


def test_ten_thousand_row_increment_is_fast_and_baseline_read_unblocked(svc):
    n = 10_000
    header = "#CHROM\tPOS\tREF\tALT\tBIG\tOTHER\n"
    rows = "\n".join(f"chr{i % 24 + 1}\t{1000 + i}\tA\tT\t0/1\t0/0" for i in range(n))
    t0 = time.perf_counter()
    v0 = svc.submit_multi_version(header + rows)
    build_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    v1 = svc.increment_multi(
        "#CHROM\tPOS\tREF\tALT\tBIG\nchr1\t1000\tA\tT\t1/1\n", v0.version_id)
    inc_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    big_rows = "\n".join(f"chr{i % 24 + 1}\t{90000 + i}\tA\tT\t0/1" for i in range(n))
    v2 = svc.increment_single(big_rows, "NEW", v1.version_id)
    big_s = time.perf_counter() - t0

    latencies = []
    for _ in range(100):
        t = time.perf_counter()
        assert svc.get_version(v0.version_id) is not None
        latencies.append(time.perf_counter() - t)

    logger.info(
        "判定依据 万行规模: 基线构建 %.3fs, 单样本修正增量 %.3fs, 万行新样本增量 %.3fs, "
        "基线读取最大 %.4fms; entries %d->%d",
        build_s, inc_s, big_s, max(latencies) * 1000,
        len(v0.result.entries), len(v2.result.entries),
    )
    assert len(v0.result.entries) == n
    assert len(v2.result.entries) == 2 * n
    assert next(e for e in v1.result.entries if e.variant.pos == 1000) \
        .genotype_by_sample["BIG"].value == "1/1"
    assert build_s < 5 and inc_s < 3 and big_s < 8
    assert max(latencies) < 0.01  # baseline reads never blocked
