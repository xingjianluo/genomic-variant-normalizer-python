"""Core record model shared by parsing, normalization and merging.

The model is intentionally pure (stdlib only, no external services) so the
whole pipeline stays usable offline.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class RejectReason(str, Enum):
    """Distinguishable reasons why an input record was rejected."""

    EMPTY_LINE = "EMPTY_LINE"
    MISSING_FIELDS = "MISSING_FIELDS"
    BAD_POSITION = "BAD_POSITION"
    BAD_ALLELE = "BAD_ALLELE"
    BAD_CHROM = "BAD_CHROM"
    BAD_GENOTYPE = "BAD_GENOTYPE"
    HEADER_SAMPLE_MISMATCH = "HEADER_SAMPLE_MISMATCH"
    DUPLICATE_SAMPLE = "DUPLICATE_SAMPLE"
    MALFORMED_LINE = "MALFORMED_LINE"


class ConflictType(str, Enum):
    """Semantic conflicts detected while merging records."""

    GENOTYPE_CONFLICT = "GENOTYPE_CONFLICT"
    REF_CONFLICT = "REF_CONFLICT"
    SAMPLE_CALL_CONFLICT = "SAMPLE_CALL_CONFLICT"


class Genotype(str, Enum):
    """Canonical genotype vocabulary shared by both input shapes."""

    MISSING = "./."      # unspecified / uncalled
    REF = "0/0"          # called reference (carries no ALT call)
    HET = "0/1"          # heterozygous ALT
    HOM = "1/1"          # homozygous ALT


@dataclass(frozen=True)
class Location:
    """Where in the submitted text a record/problem came from."""

    block: str           # input block label, e.g. "sample:S1" or "multisample"
    line_no: int         # 1-based line number within that block
    raw: str             # original raw line text


@dataclass(frozen=True)
class VariantKey:
    """Normalized identity of a variant: same key <=> equivalent writing."""

    chrom: str
    pos: int
    ref: str
    alt: str

    def render(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


@dataclass(frozen=True)
class RawRecord:
    """One parsed input line before validation/normalization."""

    chrom: str
    pos_raw: str
    ref: str
    alt: str
    genotype: str
    sample: str
    location: Location


@dataclass(frozen=True)
class VariantRecord:
    """One normalized, accepted call for one sample."""

    key: VariantKey
    sample: str
    genotype: Genotype
    location: Location


@dataclass(frozen=True)
class Rejection:
    """A rejected record with distinguishable reason and location."""

    reason: RejectReason
    location: Location
    detail: str


@dataclass(frozen=True)
class Conflict:
    """A semantic conflict between records that describe the same thing."""

    conflict_type: ConflictType
    variant: VariantKey
    sample: str | None
    detail: str
    locations: tuple[Location, ...]


@dataclass(frozen=True)
class CallsetEntry:
    """One aggregated variant across all samples."""

    variant: VariantKey
    samples: tuple[str, ...]                 # sorted, de-duplicated carriers
    genotype_by_sample: dict[str, Genotype]
    support: tuple[Location, ...]            # sorted supporting locations
    conflicts: tuple[Conflict, ...]


@dataclass(frozen=True)
class CallsetResult:
    """Outcome of processing one submission."""

    content_hash: str
    entries: list[CallsetEntry] = field(default_factory=list)
    rejected: list[Rejection] = field(default_factory=list)
    conflicts: list[Conflict] = field(default_factory=list)
    accepted_count: int = 0
