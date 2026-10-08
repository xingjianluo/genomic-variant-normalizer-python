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
    version_id: str | None = None  # set when the submission is versioned


class IncrementErrorType(str, Enum):
    """Distinguishable reasons an incremental commit cannot be applied."""

    BASELINE_NOT_FOUND = "BASELINE_NOT_FOUND"
    BASELINE_SUPERSEDED = "BASELINE_SUPERSEDED"
    EMPTY_INCREMENT = "EMPTY_INCREMENT"


@dataclass(frozen=True)
class EntryChange:
    """How one aggregated entry differs between two versions."""

    variant: VariantKey
    kind: str                              # "ADDED" | "REMOVED" | "CHANGED"
    samples_added: tuple[str, ...] = ()
    samples_removed: tuple[str, ...] = ()
    # (sample, old_genotype, new_genotype)
    genotype_changes: tuple[tuple[str, str, str], ...] = ()


@dataclass(frozen=True)
class ChangeSummary:
    """Deterministic diff of one version against its parent version."""

    added_entries: tuple[VariantKey, ...] = ()
    removed_entries: tuple[VariantKey, ...] = ()
    changed_entries: tuple[EntryChange, ...] = ()
    # sample-level carrier presence across the whole callset
    samples_added: tuple[str, ...] = ()
    samples_removed: tuple[str, ...] = ()
    new_conflicts: tuple[Conflict, ...] = ()
    resolved_conflicts: tuple[Conflict, ...] = ()
    new_rejections: tuple[Rejection, ...] = ()
    resolved_rejections: tuple[Rejection, ...] = ()
    identical: bool = False


@dataclass(frozen=True)
class CallsetSnapshot:
    """One immutable, content-addressed callset version.

    A snapshot keeps the full per-sample record state (accepted normalized
    records and rejections attributed to samples) so a later increment can
    replace one sample's calls wholesale and re-merge without the original
    text, while untouched samples stay byte-for-byte stable.
    """

    version_id: str                        # content hash of the full state
    parent_version: str | None             # version this one was built on
    root_version: str                      # lineage root (id of the root state)
    sequence: int                          # 0 for a root, >= 1 for increments
    commit_seq: int                        # monotonic commit order within lineage
    samples: tuple[str, ...]               # sorted sample names present
    touched_samples: tuple[str, ...]       # samples this commit replaced/added
    accepted: tuple[VariantRecord, ...]
    rejected: tuple[Rejection, ...]
    rejection_owners: tuple[tuple[str, ...], ...]  # parallel to rejected: owner samples
    result: CallsetResult
    changes: ChangeSummary
