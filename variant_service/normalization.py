"""Allele / coordinate normalization: equivalent writings converge to one key.

Pure functions, stdlib only. No external annotation or lift-over services are
used, so the service remains fully offline.

Equivalence rules implemented here:

* chromosome aliases collapse (``1`` == ``chr1`` == ``CHR1``);
* alleles are upper-cased and surrounding whitespace stripped;
* VCF-style left alignment: a shared leading base is removed while shifting
  POS forward, then shared trailing bases are trimmed; at least one anchor base
  is always kept (so ``A/A`` stays representable and is caught later as
  REF==ALT);
* genotypes collapse to a small canonical vocabulary.
"""
from __future__ import annotations

from .models import Genotype, VariantKey

_ALLELE_ALPHABET = frozenset("ACGT")

# Raw genotype text -> canonical genotype. Covers phased/unphased and the
# single-sample compact synonyms.
_GENOTYPE_ALIASES: dict[str, Genotype] = {}
for _token, _gt in (
    ("./.", Genotype.MISSING),
    (".|.", Genotype.MISSING),
    (".", Genotype.MISSING),
    ("", Genotype.MISSING),
    ("UNCALLED", Genotype.MISSING),
    ("0/0", Genotype.REF),
    ("0|0", Genotype.REF),
    ("REF", Genotype.REF),
    ("WT", Genotype.REF),
    ("0/1", Genotype.HET),
    ("0|1", Genotype.HET),
    ("1/0", Genotype.HET),
    ("1|0", Genotype.HET),
    ("HET", Genotype.HET),
    ("1/1", Genotype.HOM),
    ("1|1", Genotype.HOM),
    ("HOM", Genotype.HOM),
    ("HOM_ALT", Genotype.HOM),
):
    _GENOTYPE_ALIASES[_token] = _gt


def normalize_chrom(raw: str) -> str:
    """Return a canonical chromosome label, e.g. ``'chr1'`` for ``'CHR1'``.

    Raises ``ValueError`` for empty labels. Pure syntactic canonicalization;
    unknown non-empty labels (e.g. ``chrUn_KN707``) are preserved so that callers
    decide whether to reject them.
    """
    chrom = raw.strip().upper()
    if not chrom:
        raise ValueError("empty chromosome")
    if not chrom.startswith("CHR"):
        chrom = "CHR" + chrom
    return chrom


def is_valid_allele(allele: str) -> bool:
    """True for non-empty ACGTN allele strings (already upper/stripped)."""
    return bool(allele) and all(ch in _ALLELE_ALPHABET for ch in allele)


def canonicalize_genotype(raw: str) -> Genotype:
    """Map raw genotype text to the canonical enum.

    Raises ``ValueError`` for tokens that are not a recognized genotype.
    """
    token = raw.strip().upper()
    if token in _GENOTYPE_ALIASES:
        return _GENOTYPE_ALIASES[token]
    raise ValueError(f"unrecognized genotype: {raw!r}")


def normalize_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-align and trim redundant shared bases (VCF-style anchoring).

    ``ref``/``alt`` must already be stripped, upper-cased, non-empty DNA
    strings and ``pos`` a 1-based positive coordinate. Returns the shifted
    ``(pos, ref, alt)``. Each allele keeps at least one base anchor.
    """
    if pos < 1:
        raise ValueError("position must be >= 1")
    if not ref or not alt:
        raise ValueError("alleles must be non-empty")

    # 1) Strip shared suffix (never empty the alleles).
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # 2) Strip shared leading base, shifting POS right; keep one anchor base.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    return pos, ref, alt


def build_variant_key(chrom: str, pos: int, ref: str, alt: str) -> VariantKey:
    """Canonical chromosome + left-normalized alleles -> identity key."""
    canon_chrom = normalize_chrom(chrom)
    norm_pos, norm_ref, norm_alt = normalize_alleles(pos, ref, alt)
    return VariantKey(canon_chrom, norm_pos, norm_ref, norm_alt)
