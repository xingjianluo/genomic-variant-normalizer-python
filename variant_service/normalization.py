"""Allele / coordinate normalization: equivalent writings converge to one key."""
from __future__ import annotations

from .models import Genotype, VariantKey

_ALLELE_ALPHABET = frozenset("ACGT")


def normalize_chrom(raw: str) -> str:
    raise NotImplementedError


def canonicalize_genotype(raw: str) -> Genotype:
    raise NotImplementedError


def normalize_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-align and trim redundant shared bases (VCF-style anchoring)."""
    raise NotImplementedError


def build_variant_key(chrom: str, pos: int, ref: str, alt: str) -> VariantKey:
    raise NotImplementedError
