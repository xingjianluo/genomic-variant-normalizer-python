"""merge 测试：去重聚合、基因型冲突、位点冲突、稳定排序。"""
import logging

from variant_service.merge import merge_calls
from variant_service.models import Genotype, RejectReason, RejectedRecord, VariantCall

logger = logging.getLogger("test.merge")


def _call(sample, chrom, pos, ref, alt, gt, line=1):
    return VariantCall(sample_id=sample, chrom=chrom, pos=pos, ref=ref, alt=alt,
                       genotype=gt, source_line=line)


def test_duplicate_equivalent_calls_count_once():
    calls = [
        _call("S1", "chr1", 10, "A", "G", Genotype.HET),
        _call("S1", "chr1", 10, "A", "G", Genotype.HET),  # 完全重复
        _call("S2", "chr1", 10, "A", "G", Genotype.HOM_ALT),
    ]
    cs = merge_calls("id1", 1, calls, [])
    logger.info("输入摘要 3 calls(含1重复) -> variants=%d accepted=%d",
                len(cs.variants), cs.accepted_calls)
    assert len(cs.variants) == 1
    assert cs.accepted_calls == 2
    v = cs.variants[0]
    assert v.call_count == 2
    assert v.samples == ["S1", "S2"]
    assert v.conflicting is False


def test_genotype_conflict_is_marked_not_overwritten():
    calls = [
        _call("S1", "chr1", 10, "A", "G", Genotype.HET, line=2),
        _call("S1", "chr1", 10, "A", "G", Genotype.HOM_ALT, line=5),
        _call("S1", "chr1", 10, "A", "G", Genotype.UNKNOWN, line=9),  # 不参与冲突
    ]
    cs = merge_calls("id2", 1, calls, [])
    v = cs.variants[0]
    logger.info("输入摘要 同样本矛盾基因型 -> conflicting=%s types=%s summary=%s",
                v.conflicting, v.conflict_types, v.genotype_summary)
    assert v.conflicting is True
    assert "GENOTYPE_CONFLICT" in v.conflict_types
    # 两种矛盾证据都保留，没有被覆盖
    assert set(v.genotype_summary["S1"]) == {"HET", "HOM_ALT", "UNKNOWN"}


def test_hom_ref_vs_alt_is_conflict():
    calls = [
        _call("S1", "chr1", 10, "A", "G", Genotype.HOM_REF),
        _call("S1", "chr1", 10, "A", "G", Genotype.HET),
    ]
    cs = merge_calls("id3", 1, calls, [])
    assert cs.variants[0].conflicting is True


def test_site_conflict_different_ref_same_position():
    calls = [
        _call("S1", "chr1", 10, "A", "G", Genotype.HET),
        _call("S2", "chr1", 10, "AC", "AT", Genotype.HET),
        # 同 REF 不同 ALT：正常多 ALT 位点，不算位点冲突
        _call("S3", "chr1", 10, "A", "T", Genotype.HET),
    ]
    cs = merge_calls("id4", 1, calls, [])
    site_conf = [v for v in cs.variants if "SITE_CONFLICT" in v.conflict_types]
    logger.info("输入摘要 同位置不同REF + 同REF不同ALT -> 位点冲突变异数=%d", len(site_conf))
    # 该位置 REF 存在歧义，三个变异都应被明确标记（含同 REF 不同 ALT 的 A>T）
    assert len(site_conf) == 3
    multi_alt_keys = {(v.ref, v.alt) for v in cs.variants}
    assert ("A", "T") in multi_alt_keys
    at = next(v for v in cs.variants if (v.ref, v.alt) == ("A", "T"))
    assert "SITE_CONFLICT" in at.conflict_types
    assert "GENOTYPE_CONFLICT" not in at.conflict_types
    for v in site_conf:
        assert next(c for c in v.conflicts if c["type"] == "SITE_CONFLICT")["refs"] == ["A", "AC"]


def test_output_order_is_stable():
    calls = [
        _call("S1", "chrX", 1, "A", "G", Genotype.HET),
        _call("S1", "chr1", 20, "A", "G", Genotype.HET),
        _call("S1", "chr1", 10, "A", "G", Genotype.HET),
        _call("S1", "chr1", 10, "A", "T", Genotype.HET),
    ]
    cs1 = merge_calls("id", 1, calls, [])
    cs2 = merge_calls("id", 1, list(reversed(calls)), [])
    order1 = [(v.chrom, v.pos, v.ref, v.alt) for v in cs1.variants]
    order2 = [(v.chrom, v.pos, v.ref, v.alt) for v in cs2.variants]
    logger.info("判定依据: 正序/逆序输入输出顺序一致=%s, 顺序=%s", order1 == order2, order1)
    assert order1 == order2
    assert order1 == [("chr1", 10, "A", "G"), ("chr1", 10, "A", "T"),
                      ("chr1", 20, "A", "G"), ("chrX", 1, "A", "G")]
    assert cs1.model_dump() == cs2.model_dump()


def test_rejections_are_sorted_and_preserved():
    rejected = [
        RejectedRecord(source="b", line=3, reason=RejectReason.INVALID_ALLELE, detail="x", raw=""),
        RejectedRecord(source="a", line=1, reason=RejectReason.MALFORMED_LINE, detail="y", raw=""),
    ]
    cs = merge_calls("id", 0, [], rejected)
    assert [r.source for r in cs.rejected_records] == ["a", "b"]
