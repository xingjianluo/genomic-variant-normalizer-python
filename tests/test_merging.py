"""Tests for validation/normalization funnel and callset merging rules."""
import logging

from variant_service.merging import merge_callset, validate_and_normalize
from variant_service.models import (
    ConflictType,
    Genotype,
    Location,
    RawRecord,
    RejectReason,
)
from variant_service.normalization import build_variant_key

logger = logging.getLogger(__name__)


def _raw(chrom, pos, ref, alt, gt, sample, line_no=1, block="sample:S"):
    loc = Location(block=block, line_no=line_no, raw=f"{chrom} {pos} {ref} {alt} {gt}")
    return RawRecord(chrom=chrom, pos_raw=str(pos), ref=ref, alt=alt, genotype=gt,
                     sample=sample, location=loc)


def test_validation_rejects_bad_fields_with_distinct_reasons():
    raws = [
        _raw("chr1", "x", "A", "T", "0/1", "S1", 1),
        _raw("chr1", 100, "A", "U", "0/1", "S1", 2),
        _raw("chr1", 100, "A", "A", "0/1", "S1", 3),
        _raw("", 100, "A", "T", "0/1", "S1", 4),
        _raw("chr1", 100, "A", "T", "9/9", "S1", 5),
    ]
    accepted, rejected = validate_and_normalize(raws)
    reasons = [r.reason for r in rejected]
    logger.info("输入摘要 5 条畸形记录 -> 原因=%s 行号=%s",
                [r.reason.value for r in rejected], [r.location.line_no for r in rejected])
    assert accepted == []
    assert reasons == [
        RejectReason.BAD_POSITION,
        RejectReason.BAD_ALLELE,
        RejectReason.BAD_ALLELE,
        RejectReason.BAD_CHROM,
        RejectReason.BAD_GENOTYPE,
    ]
    # location preserved per record
    assert [r.location.line_no for r in rejected] == [1, 2, 3, 4, 5]


def test_equivalent_writings_aggregate_into_one_entry():
    raws = [
        _raw("1", 100, "ACG", "AT", "0/1", "S1", 1),
        _raw("chr1", 101, "CG", "T", "1/1", "S2", 2),
    ]
    accepted, rejected = validate_and_normalize(raws)
    assert rejected == [] and len(accepted) == 2
    result = merge_callset("h", accepted, [])
    logger.info("输入摘要 等价写法两条 -> entries=%d carriers=%s",
                len(result.entries), [e.samples for e in result.entries])
    assert len(result.entries) == 1
    entry = result.entries[0]
    assert entry.samples == ("S1", "S2")
    assert entry.genotype_by_sample["S1"] is Genotype.HET
    assert entry.genotype_by_sample["S2"] is Genotype.HOM


def test_duplicate_identical_calls_do_not_double_count():
    raws = [
        _raw("chr1", 100, "A", "T", "0/1", "S1", 1),
        _raw("chr1", 100, "A", "T", "0/1", "S1", 2),
        _raw("1", 100, "A", "T", "0/1", "S1", 3),
    ]
    accepted, _ = validate_and_normalize(raws)
    result = merge_callset("h", accepted, [])
    entry = result.entries[0]
    assert entry.samples == ("S1",)
    assert len(entry.support) == 3  # evidence kept, sample counted once
    assert result.conflicts == []


def test_genotype_conflict_is_marked_not_resolved():
    raws = [
        _raw("chr1", 100, "A", "T", "0/1", "S1", 1),
        _raw("chr1", 100, "A", "T", "1/1", "S1", 2),
    ]
    accepted, _ = validate_and_normalize(raws)
    result = merge_callset("h", accepted, [])
    logger.info("输入摘要 同样本杂合/纯合矛盾 -> conflicts=%s",
                [c.conflict_type.value for c in result.conflicts])
    assert len(result.entries) == 1
    assert [c.conflict_type for c in result.conflicts] == [ConflictType.GENOTYPE_CONFLICT]
    assert result.entries[0].conflicts[0].sample == "S1"


def test_ref_conflict_between_hom_alt_and_hom_ref():
    raws = [
        _raw("chr1", 100, "A", "T", "1/1", "S1", 1),
        _raw("chr1", 100, "A", "T", "0/0", "S1", 2),
    ]
    accepted, _ = validate_and_normalize(raws)
    result = merge_callset("h", accepted, [])
    types = {c.conflict_type for c in result.conflicts}
    logger.info("输入摘要 纯合ALT/纯合REF矛盾 -> %s", [t.value for t in types])
    assert ConflictType.REF_CONFLICT in types
    # 0/0 is not a carrier
    assert result.entries[0].samples == ("S1",)
    assert result.entries[0].genotype_by_sample["S1"] is Genotype.HOM


def test_distinct_alts_same_sample_same_site_flagged_on_both_entries():
    raws = [
        _raw("chr1", 100, "A", "T", "0/1", "S1", 1),
        _raw("chr1", 100, "A", "G", "0/1", "S1", 2),
    ]
    accepted, _ = validate_and_normalize(raws)
    result = merge_callset("h", accepted, [])
    assert len(result.entries) == 2  # neither ALT dropped
    flagged = [e for e in result.entries if e.conflicts]
    assert len(flagged) == 2
    assert all(c.conflict_type is ConflictType.SAMPLE_CALL_CONFLICT for c in result.conflicts)
    # different samples with different ALTs is NOT a conflict
    raws2 = [
        _raw("chr1", 100, "A", "T", "1/1", "S1", 1),
        _raw("chr1", 100, "A", "G", "1/1", "S2", 2),
    ]
    accepted2, _ = validate_and_normalize(raws2)
    result2 = merge_callset("h2", accepted2, [])
    logger.info("判定依据: 跨样本不同 ALT conflicts=%s", len(result2.conflicts))
    assert result2.conflicts == []
    assert {e.samples for e in result2.entries} == {("S1",), ("S2",)}


def test_output_order_is_deterministic():
    raws = [
        _raw("chr10", 5, "A", "T", "0/1", "S1", 1),
        _raw("chr2", 5, "A", "T", "0/1", "S1", 2),
        _raw("chr1", 5, "A", "T", "0/1", "S1", 3),
        _raw("chr1", 2, "A", "T", "0/1", "S1", 4),
        _raw("chrX", 1, "A", "T", "0/1", "S1", 5),
    ]
    import random
    results = []
    for seed in range(5):
        shuffled = raws[:]
        random.Random(seed).shuffle(shuffled)
        accepted, _ = validate_and_normalize(shuffled)
        r = merge_callset("h", accepted, [])
        results.append([(e.variant.chrom, e.variant.pos) for e in r.entries])
    logger.info("输入摘要 乱序 5 次 -> 首次顺序=%s", results[0])
    assert all(r == results[0] for r in results)
    assert results[0] == [("CHR1", 2), ("CHR1", 5), ("CHR2", 5), ("CHR10", 5), ("CHRX", 1)]


def test_ref_only_records_are_not_carriers_but_kept():
    raws = [_raw("chr1", 100, "A", "T", "0/0", "S1", 1)]
    accepted, _ = validate_and_normalize(raws)
    result = merge_callset("h", accepted, [])
    assert result.entries[0].samples == ()
    assert result.entries[0].genotype_by_sample["S1"] is Genotype.REF
