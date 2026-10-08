"""Incremental update tests: add/correct samples, versioning, change summary.

Covers the service layer: root creation, per-sample wholesale replacement
(never evidence stacking), deterministic content-addressed versions,
idempotent replay, the change-summary taxonomy, structural rejection
attribution, and optimistic-concurrency behaviour for disjoint and
overlapping concurrent increments.
"""
import logging
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from variant_service.models import ConflictType, RejectReason
from variant_service.service import BaselineError, CallsetService

logger = logging.getLogger(__name__)


def _entry(snap, render):
    return {e.variant.render(): e for e in snap.result.entries}[render]


def test_root_version_is_stable_and_serialized_state(service):
    root = service.submit_root_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
    )
    again = service.submit_root_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
    )
    logger.info(
        "判定依据[根版本幂等]: 相同文本两次建根 -> 同一对象 %s, version=%s 前 8 位, samples=%s",
        root is again, root.version_id[:8], root.samples,
    )
    assert root is again  # content-addressed, object reused
    assert root.sequence == 0
    assert root.parent_version is None
    assert root.root_version == root.version_id
    assert service.get_version(root.version_id) is root


def test_increment_adds_new_sample_without_touching_existing(service):
    root = service.submit_root_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\nchr1\t100\tA\tT\t0/1\t0/1\n"
    )
    v2 = service.submit_increment_single("chr2 200 G C 0/1\n", "S3", root.version_id)
    logger.info(
        "判定依据[新增样本]: untouched S1/S2 条目不变=%s, 新条目=%s, samples_added=%s",
        [e.samples for e in v2.result.entries if e.variant.render() == "CHR1:100:A:T"],
        [v.render() for v in v2.changes.added_entries], v2.changes.samples_added,
    )
    # old entry keeps exact carriers and order; nothing removed/reordered
    assert [e.variant.render() for e in v2.result.entries] == [
        "CHR1:100:A:T", "CHR2:200:G:C",
    ]
    assert _entry(v2, "CHR1:100:A:T").samples == ("S1", "S2")
    assert v2.changes.added_entries[0].render() == "CHR2:200:G:C"
    assert v2.changes.samples_added == ("S3",)
    assert v2.changes.removed_entries == ()
    assert v2.parent_version == root.version_id
    assert v2.root_version == root.root_version
    assert v2.sequence == 1


def test_correction_replaces_sample_not_stacks_evidence(service):
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    # correction: same sample now says HOM; must replace, not stack HET+HOM
    v2 = service.submit_increment_single("chr1 100 A T 1/1\n", "S1", root.version_id)
    entry = _entry(v2, "CHR1:100:A:T")
    logger.info(
        "判定依据[修正非叠加]: 修正后基因型=%s, GENOTYPE_CONFLICT 数=%d (替换而非新旧证据叠加)",
        entry.genotype_by_sample["S1"].value,
        len([c for c in v2.result.conflicts
             if c.conflict_type is ConflictType.GENOTYPE_CONFLICT]),
    )
    assert entry.genotype_by_sample["S1"].value == "1/1"
    assert all(c.conflict_type is not ConflictType.GENOTYPE_CONFLICT
               for c in v2.result.conflicts)
    change = v2.changes.changed_entries[0]
    assert change.kind == "CHANGED"
    assert change.genotype_changes == (("S1", "0/1", "1/1"),)


def test_correction_can_remove_entries_and_carriers(service):
    root = service.submit_root_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\t0/1\n"
        "chr2\t200\tG\tC\t0/1\t./.\n"
    )
    # S1 retracts both calls; its carrier presence disappears, S2 stays
    v2 = service.submit_increment_single("", "S1", root.version_id)
    logger.info(
        "判定依据[样本撤回归零]: removed_entries=%s, 残留样本=%s",
        [v.render() for v in v2.changes.removed_entries],
        [e.samples for e in v2.result.entries],
    )
    assert [v.render() for v in v2.changes.removed_entries] == ["CHR2:200:G:C"]
    left = _entry(v2, "CHR1:100:A:T")
    assert left.samples == ("S2",)
    assert "S1" not in left.genotype_by_sample
    # S1 was an ALT carrier before, so it leaves the global carrier set
    assert v2.changes.samples_removed == ("S1",)


def test_change_summary_new_and_resolved_conflicts(service):
    root = service.submit_root_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\nchr1\t100\tA\tT\t1/1\t0/1\n"
    )
    # S1 introduces a second ALT at the same site -> SAMPLE_CALL_CONFLICT
    v2 = service.submit_increment_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\n"
        "chr1\t100\tA\tT\t1/1\n"
        "chr1\t100\tA\tG\t0/1\n",
        root.version_id,
    )
    new_types = {c.conflict_type for c in v2.changes.new_conflicts}
    logger.info("判定依据[新冲突]: new_conflicts=%s", {t.value for t in new_types})
    assert ConflictType.SAMPLE_CALL_CONFLICT in new_types

    # S1 corrects back to one ALT -> conflict resolves
    v3 = service.submit_increment_single("chr1 100 A T 1/1\n", "S1", v2.version_id)
    resolved_types = {c.conflict_type for c in v3.changes.resolved_conflicts}
    logger.info(
        "判定依据[冲突消除]: resolved=%s, 现存冲突=%d",
        {t.value for t in resolved_types}, len(v3.result.conflicts),
    )
    assert ConflictType.SAMPLE_CALL_CONFLICT in resolved_types
    assert v3.result.conflicts == []


def test_conflict_with_other_sample_still_marked_on_correction(service):
    # S2 calls HOM at a site where it also (other ALT via S1's site?) Build:
    # S1: A>T HET; correction makes S1 A>T HOM *and* A>G HET (ambiguous site).
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    v2 = service.submit_increment_multi(
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t1/1\t0/1\n"
        "chr1\t100\tA\tG\t1/1\t./.\n",
        root.version_id,
    )
    types = {c.conflict_type for c in v2.result.conflicts}
    logger.info("判定依据[增量沿用冲突规则]: 冲突类型=%s", {t.value for t in types})
    assert ConflictType.SAMPLE_CALL_CONFLICT in types


def test_malformed_increment_keeps_distinguishable_rejection(service):
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    v2 = service.submit_increment_single(
        "chr1 x A T 0/1\nbad row\nchr1 100 A U 0/1\n", "S2", root.version_id
    )
    reasons = [(r.reason, r.location.line_no) for r in v2.result.rejected]
    logger.info("判定依据[畸形增量不入账]: rejected=%s (旧 S1 条目保留=%d)",
                reasons,
                len([e for e in v2.result.entries if e.variant.render() == "CHR1:100:A:T"]))
    assert reasons == [
        (RejectReason.BAD_POSITION, 1),
        (RejectReason.MISSING_FIELDS, 2),
        (RejectReason.BAD_ALLELE, 3),
    ]
    # bad sample added nothing, old entry intact
    assert _entry(v2, "CHR1:100:A:T").samples == ("S1",)

    # correcting S2 with valid text clears its prior rejections
    v3 = service.submit_increment_single("chr3 300 C G 1/1\n", "S2", v2.version_id)
    assert [r.reason.value for r in v3.changes.resolved_rejections].count(
        "BAD_POSITION") == 1
    assert v3.result.rejected == []


def test_column_mismatch_in_increment_is_rejected_not_swallowed(service):
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    v2 = service.submit_increment_multi(
        "#CHROM\tPOS\tREF\tALT\tS2\tS3\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
        "chr2\t200\tG\tC\t0/1\n",  # missing one sample column
        root.version_id,
    )
    reasons = {r.reason for r in v2.result.rejected}
    logger.info("判定依据[样本列错位]: 拒绝类型=%s", {r.value for r in reasons})
    assert RejectReason.HEADER_SAMPLE_MISMATCH in reasons


def test_same_baseline_same_increment_is_reproducible(service):
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    payload = ("#CHROM\tPOS\tREF\tALT\tS2\n"
               "chr1\t100\tA\tT\t1/1\nchr3\t300\tC\tG\t0/1\n")
    a = service.submit_increment_multi(payload, root.version_id)
    b = service.submit_increment_multi(payload, root.version_id)
    c = service.submit_increment_multi(payload, root.version_id)
    logger.info(
        "判定依据[可复现]: 三次提交 version 一致=%s, 摘要 added=%s 一致",
        a.version_id == b.version_id == c.version_id,
        [v.render() for v in a.changes.added_entries],
    )
    assert a.version_id == b.version_id == c.version_id
    # same baseline + same increment => byte-identical summary every time
    assert a.changes == b.changes == c.changes
    assert b.result.entries == c.result.entries == a.result.entries
    assert b.result.rejected == c.result.rejected == a.result.rejected
    assert b.result.conflicts == c.result.conflicts == a.result.conflicts


def test_resubmitting_current_state_is_noop(service):
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    v2 = service.submit_increment_single("chr1 100 A T 1/1\n", "S1", root.version_id)
    v3 = service.submit_increment_single("chr1 100 A T 1/1\n", "S1", v2.version_id)
    logger.info(
        "判定依据[幂等无操作]: version 不变=%s, identical=%s, sequence 停留=%d",
        v3.version_id == v2.version_id, v3.changes.identical, v3.sequence,
    )
    assert v3.version_id == v2.version_id
    assert v3.changes.identical is True
    assert v3.sequence == v2.sequence


def test_unknown_and_stale_baseline_have_distinguishable_errors(service):
    with pytest.raises(BaselineError) as ei:
        service.submit_increment_single("chr1 100 A T 0/1\n", "S9", "deadbeef")
    assert ei.value.error_type.value == "BASELINE_NOT_FOUND"
    logger.info("判定依据[未知基线]: reason=%s", ei.value.error_type.value)

    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    # first valid correction advances the lineage head
    service.submit_increment_single("chr1 100 A T 1/1\n", "S1", root.version_id)
    with pytest.raises(BaselineError) as ei2:
        # overlapping correction against the stale root cannot be merged
        service.submit_increment_single("chr1 100 A G 0/1\n", "S1", root.version_id)
    err = ei2.value
    logger.info(
        "判定依据[过期基线+样本重叠]: reason=%s, current_version 存在=%s, conflicting=%s",
        err.error_type.value, bool(err.current_version), err.conflicting_samples,
    )
    assert err.error_type.value == "BASELINE_SUPERSEDED"
    assert err.conflicting_samples == ("S1",)
    assert err.current_version

    # disjoint-sample update against the stale root is auto-rebased onto head
    rebased = service.submit_increment_single(
        "chr9 900 T A 0/1\n", "S2", root.version_id
    )
    assert "S1" in rebased.samples and "S2" in rebased.samples
    logger.info("判定依据[不相交样本自动重放]: 重放后 samples=%s 全部保留", rebased.samples)


def test_concurrent_disjoint_increments_lose_no_samples():
    svc = CallsetService()
    root = svc.submit_root_single("chr1 100 A T 0/1\n", "S0")

    def job(i):
        return svc.submit_increment_single(
            f"chr{i % 24 + 1}\t{2000 + i}\tG\tC\t0/1\n", f"S{i}", root.version_id
        )

    with ThreadPoolExecutor(max_workers=16) as pool:
        snaps = list(pool.map(job, range(1, 65)))
    final = svc.get_version(svc._lineage[root.root_version]["head"])
    expected = {f"S{i}" for i in range(65)}
    logger.info(
        "判定依据[并发不丢样本]: 64 线程并发后 head samples=%d, 与期望一致=%s, entries=%d",
        len(final.samples), set(final.samples) == expected, len(final.result.entries),
    )
    assert set(final.samples) == expected
    assert len(final.result.entries) == 65

    # determinism: sequential same payloads reach the identical content version
    seq = CallsetService()
    r0 = seq.submit_root_single("chr1 100 A T 0/1\n", "S0")
    cur = r0
    for i in range(1, 65):
        cur = seq.submit_increment_single(
            f"chr{i % 24 + 1}\t{2000 + i}\tG\tC\t0/1\n", f"S{i}", r0.version_id
        )
    assert cur.version_id == final.version_id


def test_concurrent_overlapping_increments_are_rejected_distinguishably():
    svc = CallsetService()
    root = svc.submit_root_single("chr1 100 A T 0/1\n", "S1")

    def job(i):
        try:
            svc.submit_increment_single(
                f"chr1 100 A T {'1/1' if i % 2 else '0/1'}\n",
                "S1", root.version_id,
            )
            return "ok"
        except BaselineError as err:
            return err.error_type.value

    with ThreadPoolExecutor(max_workers=32) as pool:
        outcomes = list(pool.map(job, range(48)))
    rejected = [o for o in outcomes if o != "ok"]
    logger.info(
        "判定依据[并发重叠可复现拒绝]: ok=%d, 拒绝=%d, 拒绝原因集合=%s",
        outcomes.count("ok"), len(rejected), set(rejected),
    )
    assert set(rejected) <= {"BASELINE_SUPERSEDED"}
    assert outcomes.count("ok") >= 1
    # head is always a well-formed callset (never corrupted/interleaved)
    head = svc.get_version(svc._lineage[root.root_version]["head"])
    assert head.samples == ("S1",)
    assert len(head.result.entries) == 1


def test_increment_accepts_vcf_shaped_multi_payload(service):
    root = service.submit_root_multi(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n"
    )
    v2 = service.submit_increment_multi(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n"
        "chr2\t200\t.\tG\tC\t.\t.\t.\tGT:DP\t1/1:30\n",
        root.version_id,
    )
    by_key = {e.variant.render(): e for e in v2.result.entries}
    logger.info(
        "判定依据[增量兼容 VCF 形态]: 条目=%s, 拒绝=%d, S2 携带=%s",
        sorted(by_key), len(v2.result.rejected), by_key["CHR2:200:G:C"].samples,
    )
    assert not v2.result.rejected
    assert by_key["CHR1:100:A:T"].samples == ("S1",)
    assert by_key["CHR2:200:G:C"].samples == ("S2",)
    assert by_key["CHR2:200:G:C"].genotype_by_sample["S2"].value == "1/1"


def test_stale_baseline_distinct_corrections_are_rejected_not_merged(service):
    """Stale baseline + overlapping sample + genuinely new target state => 409.

    After the lineage advances through several *other* changes, a late commit
    against the old root cannot be idempotent (its target state does not exist)
    and cannot be auto-merged (it touches an already-changed sample), so it must
    be deterministically rejected rather than silently overwriting anyone.
    """
    root = service.submit_root_single("chr1 100 A T 0/1\n", "S1")
    v1 = service.submit_increment_single("chr1 100 A T 1/1\n", "S1", root.version_id)
    service.submit_increment_single("chr2 200 G C 0/1\n", "S2", v1.version_id)
    # A unique correction against the stale root: target state never seen.
    with pytest.raises(BaselineError) as ei:
        service.submit_increment_single("chr1 100 A G 0/1\n", "S1", root.version_id)
    err = ei.value
    logger.info(
        "判定依据[陈旧基线新结论拒绝]: reason=%s, conflicting=%s, current 非空=%s",
        err.error_type.value, err.conflicting_samples, bool(err.current_version),
    )
    assert err.error_type.value == "BASELINE_SUPERSEDED"
    assert err.conflicting_samples == ("S1",)
    assert err.current_version
    # retrying against the current head succeeds (deterministic recovery path)
    recovered = service.submit_increment_single(
        "chr1 100 A G 0/1\n", "S1", err.current_version
    )
    assert "S1" in recovered.samples


def test_ten_thousand_record_increment_is_fast_and_reads_unblocked():
    svc = CallsetService()
    root_rows = "\n".join(
        f"chr{i % 24 + 1}\t{1000 + i}\tA\tT\t0/1" for i in range(10_000)
    )
    root = svc.submit_root_multi("#CHROM\tPOS\tREF\tALT\tS0\n" + root_rows)

    barrier = []
    import threading
    ev = threading.Event()

    def big_correct():
        ev.wait()
        text = "#CHROM\tPOS\tREF\tALT\tS0\n" + "\n".join(
            f"chr{i % 24 + 1}\t{1000 + i}\tA\tT\t1/1" for i in range(10_000)
        )
        return svc.submit_increment_multi(text, root.version_id)

    t = threading.Thread(target=lambda: barrier.append(big_correct()))
    t.start()
    ev.set()
    start = time.perf_counter()
    # version reads stay responsive while the heavy rebuild runs
    got = svc.get_version(root.version_id)
    read_elapsed = time.perf_counter() - start
    t.join()
    updated = barrier[0]
    total = time.perf_counter() - start
    logger.info(
        "判定依据[万行性能+读不阻塞]: 重建后 entries=%d, 读旧版本耗时 %.4fs, 总耗时 %.3fs",
        len(updated.result.entries), read_elapsed, total,
    )
    assert got is root
    assert read_elapsed < 0.1
    assert len(updated.result.entries) == 10_000
    assert total < 10.0
