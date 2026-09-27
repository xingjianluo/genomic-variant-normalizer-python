"""Tests for structural parsing of both input shapes."""
import logging

from variant_service.models import RejectReason
from variant_service.parsing import parse_multi_sample, parse_single_sample

logger = logging.getLogger(__name__)


def test_single_sample_basic():
    text = "# comment\nchr1 100 A T 0/1\nchr2 200 ACG AT hom\n"
    raws, rejected = parse_single_sample(text, "S1")
    logger.info("输入摘要: %r -> %d 条原始记录, %d 拒绝", text, len(raws), len(rejected))
    assert len(raws) == 2
    assert rejected == []
    assert raws[0].sample == "S1"
    assert raws[0].genotype == "0/1"
    assert raws[0].location.line_no == 2
    assert raws[1].genotype == "hom"


def test_single_sample_missing_fields_not_silently_dropped():
    text = "chr1 100 A\nchr1 100 A T"
    raws, rejected = parse_single_sample(text, "S1")
    logger.info("输入摘要: %r -> rejected=%s", text, [(r.reason.value, r.location.line_no) for r in rejected])
    assert raws[0].pos_raw == "100"
    assert len(rejected) == 1
    assert rejected[0].reason is RejectReason.MISSING_FIELDS
    assert rejected[0].location.line_no == 1
    assert rejected[0].location.block == "sample:S1"


def test_single_sample_too_many_fields():
    raws, rejected = parse_single_sample("chr1 100 A T 0/1 EXTRA", "S1")
    assert raws == []
    assert rejected[0].reason is RejectReason.MALFORMED_LINE


def test_single_sample_no_genotype_is_missing():
    raws, rejected = parse_single_sample("chr1 100 A T", "S1")
    assert rejected == []
    assert raws[0].genotype == ""


def test_multi_sample_compact():
    text = (
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
        "chr2\t200\tG\tC\t.\t0/0\n"
    )
    raws, rejected = parse_multi_sample(text)
    logger.info("输入摘要 compact 多样本 -> %d 记录 %d 拒绝", len(raws), len(rejected))
    assert len(raws) == 4
    assert rejected == []
    assert {(r.sample, r.genotype) for r in raws} >= {("S1", "0/1"), ("S2", "1/1"), ("S2", "0/0")}


def test_multi_sample_vcf_format():
    text = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP\t0/1:12\t1/1:30\n"
    )
    raws, rejected = parse_multi_sample(text)
    assert len(raws) == 2
    assert rejected == []
    assert {r.genotype for r in raws} == {"0/1", "1/1"}


def test_multi_sample_column_count_mismatch():
    text = (
        "#CHROM\tPOS\tREF\tALT\tS1\tS2\n"
        "chr1\t100\tA\tT\t0/1\n"              # one sample cell short
        "chr2\t200\tG\tC\t0/0\t0/1\n"
    )
    raws, rejected = parse_multi_sample(text)
    logger.info("输入摘要 列数错位 -> reasons=%s", [r.reason.value for r in rejected])
    assert len(raws) == 2  # good row still parsed, bad row not silently dropped
    assert len(rejected) == 1
    assert rejected[0].reason is RejectReason.HEADER_SAMPLE_MISMATCH
    assert rejected[0].location.line_no == 2


def test_multi_sample_vcf_missing_gt_field():
    text = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tT\t.\t.\t.\tDP\t12\n"
    )
    raws, rejected = parse_multi_sample(text)
    assert raws == []
    assert rejected[0].reason is RejectReason.MALFORMED_LINE


def test_multi_sample_duplicate_sample_names():
    text = (
        "#CHROM\tPOS\tREF\tALT\tS1\tS1\n"
        "chr1\t100\tA\tT\t0/1\t1/1\n"
    )
    raws, rejected = parse_multi_sample(text)
    assert [r.reason for r in rejected] == [RejectReason.DUPLICATE_SAMPLE]
    names = sorted({r.sample for r in raws})
    logger.info("输入摘要 重名样本 -> 样本名=%s", names)
    assert names == ["S1", "S1_2"]


def test_multi_sample_without_header_infers_names():
    text = "chr1\t100\tA\tT\t0/1\t1/1\nchr2\t200\tG\tC\t0/0\t.\n"
    raws, rejected = parse_multi_sample(text)
    assert rejected == []
    assert {r.sample for r in raws} == {"sample_1", "sample_2"}
    # a later row with different count is a mismatch, not a silent resize
    text_bad = text + "chr3\t300\tA\tG\t0/1\n"
    raws2, rejected2 = parse_multi_sample(text_bad)
    assert rejected2[0].reason is RejectReason.HEADER_SAMPLE_MISMATCH
    assert len(raws2) == 4


def test_blank_and_comment_lines_are_skipped():
    raws, rejected = parse_single_sample("\n   \n# hi\nchr1 100 A T", "S1")
    assert len(raws) == 1 and rejected == []
