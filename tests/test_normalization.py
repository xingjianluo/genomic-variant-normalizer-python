"""Tests for chromosome / genotype / allele normalization."""
import logging

import pytest

from variant_service.models import Genotype
from variant_service.normalization import (
    build_variant_key,
    canonicalize_genotype,
    is_valid_allele,
    normalize_alleles,
    normalize_chrom,
)

logger = logging.getLogger(__name__)


@pytest.mark.parametrize(
    "raw,expected",
    [("1", "CHR1"), ("chr1", "CHR1"), ("CHR1", "CHR1"), (" X ", "CHRX"), ("chrM", "CHRM")],
)
def test_chrom_aliases_converge(raw, expected):
    got = normalize_chrom(raw)
    logger.info("输入摘要 chrom=%r -> %s (期望 %s)", raw, got, expected)
    assert got == expected


def test_empty_chrom_rejected():
    with pytest.raises(ValueError):
        normalize_chrom("  ")


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("0/1", Genotype.HET),
        ("1|0", Genotype.HET),
        ("1/1", Genotype.HOM),
        ("hom", Genotype.HOM),
        ("0/0", Genotype.REF),
        (".", Genotype.MISSING),
        ("", Genotype.MISSING),
    ],
)
def test_genotype_canonicalization(raw, expected):
    got = canonicalize_genotype(raw)
    logger.info("输入摘要 gt=%r -> %s", raw, got.value)
    assert got is expected


def test_bad_genotype_raises():
    with pytest.raises(ValueError):
        canonicalize_genotype("2/3")


def test_allele_validation():
    assert is_valid_allele("ACGT")
    assert not is_valid_allele("")
    assert not is_valid_allele("ACGU")


@pytest.mark.parametrize(
    "pos,ref,alt,epos,eref,ealt",
    [
        # SNV unchanged
        (100, "A", "T", 100, "A", "T"),
        # shared leading base shifts POS forward (classic VCF anchoring)
        (100, "ACG", "AT", 101, "CG", "T"),
        (101, "CG", "T", 101, "CG", "T"),
        # shared leading base shifts POS forward
        (100, "AC", "AT", 101, "C", "T"),
        # deletion: keep one anchor
        (5, "ATC", "A", 5, "ATC", "A"),
    ],
)
def test_allele_normalization(pos, ref, alt, epos, eref, ealt):
    got = normalize_alleles(pos, ref, alt)
    logger.info("输入摘要 %d %s>%s -> %s", pos, ref, alt, got)
    assert got == (epos, eref, ealt)


def test_equivalent_writings_build_same_key():
    # Two writings of the same indel with redundant anchors must converge.
    k1 = build_variant_key("1", 100, "ACG", "AT")
    k2 = build_variant_key("chr1", 101, "CG", "T")
    logger.info("等价写法判定依据: key1=%s key2=%s", k1.render(), k2.render())
    assert k1 == k2


def test_distinct_variants_build_different_keys():
    assert build_variant_key("1", 100, "A", "T") != build_variant_key("1", 100, "A", "G")
    assert build_variant_key("1", 100, "A", "T") != build_variant_key("2", 100, "A", "T")


def test_key_is_frozen_and_hashable():
    k = build_variant_key("chr1", 10, "A", "T")
    assert hash(k) == hash(build_variant_key("1", 10, "A", "T"))
    with pytest.raises(Exception):
        k.pos = 11  # type: ignore[misc]
