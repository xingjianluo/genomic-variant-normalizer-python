"""normalization 测试：等价写法归一 + 畸形输入分类。"""
import logging

from variant_service.models import Genotype, RejectReason
from variant_service.normalization import (
    normalize_allele,
    normalize_chrom,
    parse_input,
    parse_single,
    parse_vcf,
)

logger = logging.getLogger("test.normalization")


def test_chrom_equivalent_spellings():
    cases = {"1": "chr1", "CHR1": "chr1", "chr13": "chr13", "X": "chrX", "chrm": "chrMT"}
    for raw, expected in cases.items():
        got = normalize_chrom(raw)
        logger.info("输入摘要 chrom=%r -> 判定 %s (期望 %s)", raw, got, expected)
        assert got == expected
    assert normalize_chrom("chrZ") is None


def test_allele_trim_rules():
    # 先去前缀再去后缀（VCF 左移）: ACG / ATG@5 -> C / T@6
    assert normalize_allele("ACG", "ATG", 5) == (6, "C", "T")
    # 公共前缀 trim（带左移）: CA / CTA@100 -> A / TA@101
    got = normalize_allele("CA", "CTA", 100)
    logger.info("输入摘要 CA>CTA@100 => %s", got)
    assert got == (101, "A", "TA")
    # 大小写归一
    assert normalize_allele("a", "g", 7) == (7, "A", "G")
    # 非法
    assert normalize_allele("A", "A", 1) is None
    assert normalize_allele("A", "<DEL>", 1) is None
    assert normalize_allele("", "G", 1) is None


def test_single_equivalent_records_converge():
    # 同一位置的等价写法（大小写 / 染色体前缀 / 冗余锚定）必须收敛
    text = "\n".join([
        "S1\t1\t100\tAC\tGC\tHET",    # trim 公共后缀 -> A/G@100
        "S1\tCHR1\t100\ta\tg\t0/1",
        "S1\tchr1\t100\tCA\tCTA\t1/1",  # trim 公共前缀 -> A/TA@101，不同变异
    ])
    result = parse_single("", text)
    logger.info("输入摘要 3 行文本，accept=%d reject=%d", len(result.calls), len(result.rejected))
    assert len(result.rejected) == 0
    keys = {(c.chrom, c.pos, c.ref, c.alt) for c in result.calls}
    logger.info("判定归一键集合=%s", keys)
    assert keys == {("chr1", 100, "A", "G"), ("chr1", 101, "A", "TA")}


def test_malformed_single_lines_are_classified():
    text = "\n".join([
        "S1\tchr1\t100\tA",          # 列数不足
        "S1\tchr1\tx\tA\tG",         # 位置非整数
        "S1\tchrZ\t100\tA\tG",       # 染色体无法识别
        "S1\tchr1\t100\tA\t<DEL>",   # 非法等位
        "S1\tchr1\t100\tA\tG\tFOO/BAR",  # 无法识别的基因型
        "S1\tchr1\t0\tA\tG",         # 位置越界
        "",
        "# 注释行",
    ])
    result = parse_single("", text)
    reasons = [r.reason for r in result.rejected]
    logger.info("输入摘要 畸形行 6 条，拒绝原因=%s", reasons)
    assert RejectReason.WRONG_COLUMN_COUNT in reasons
    assert RejectReason.INVALID_POSITION in reasons
    assert RejectReason.INVALID_CHROM in reasons
    assert RejectReason.INVALID_ALLELE in reasons
    assert RejectReason.INVALID_GENOTYPE in reasons
    assert all(r.line >= 1 for r in result.rejected)
    logger.info("判定依据: 每条拒绝带行号与原因，互不混淆")


def test_vcf_multisample_expansion_and_mismatch():
    vcf = "\n".join([
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1",
        "chr1\t20\t.\tC\tT,A\t.\t.\t.\tGT\t0/2\t0/0",
        "chr2\t30\t.\tG\tA\t.\t.\t.\tGT\t0/1",  # 样本列对不上
    ])
    result = parse_vcf(vcf)
    logger.info("输入摘要 VCF 3 数据行 calls=%d rejected=%d", len(result.calls), len(result.rejected))
    # 行3: S1 HET + S2 HOM_ALT = 2；行4: S1 展开 index2(=A) 1 条, S2 HOM_REF 对2个ALT 2 条 = 3
    assert len(result.calls) == 5
    s1_calls = [c for c in result.calls if c.sample_id == "S1"]
    assert {c.genotype for c in s1_calls} == {Genotype.HET}
    assert {(c.pos, c.alt) for c in s1_calls} == {(10, "G"), (20, "A")}
    mismatch = [r for r in result.rejected if r.reason == RejectReason.SAMPLE_COLUMN_MISMATCH]
    assert len(mismatch) == 1 and mismatch[0].line == 5
    logger.info("判定依据: 行5 样本值列数 1 != 表头 2 -> SAMPLE_COLUMN_MISMATCH")


def test_vcf_sites_only_and_no_call():
    vcf = "\n".join([
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t10\t.\tA\tG\t.\t.\t.",
    ])
    r = parse_vcf(vcf)
    assert len(r.calls) == 1 and r.calls[0].sample_id == "UNKNOWN"
    vcf2 = "\n".join([
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t./.",
    ])
    r2 = parse_vcf(vcf2)
    assert r2.calls[0].genotype == Genotype.UNKNOWN
    logger.info("判定: 无样本列->UNKNOWN 样本; ./. -> UNKNOWN 基因型")


def test_vcf_missing_header_and_duplicate_sample():
    bad = "chr1\t10\t.\tA\tG\t.\t.\t.\n"
    r1 = parse_input("x", "vcf", bad)
    assert r1.rejected[0].reason == RejectReason.MALFORMED_LINE
    dup = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1\n"
    r2 = parse_vcf(dup)
    assert r2.rejected[0].reason == RejectReason.DUPLICATE_SAMPLE_HEADER
    logger.info("判定: 缺表头=MALFORMED_LINE; 重名样本=DUPLICATE_SAMPLE_HEADER")


def test_empty_input_rejected():
    r = parse_input("S1", "single", "   \n")
    assert r.rejected[0].reason == RejectReason.EMPTY_INPUT
