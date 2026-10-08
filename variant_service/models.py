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
    """A rejected record with distinguishable reason and location.

    ``sample`` is set when the rejection can be attributed to one concrete
    sample cell (e.g. an unrecognized genotype in one multi-sample column);
    row-level structural problems shared by a whole submission leave it None.
    """

    reason: RejectReason
    location: Location
    detail: str
    sample: str | None = None


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


class EntryChangeType(str, Enum):
    """How an entry changed between two callset versions."""

    ADDED = "ADDED"          # variant entry did not exist in the baseline
    REMOVED = "REMOVED"      # present in baseline, gone after the increment
    MODIFIED = "MODIFIED"    # same variant, carriers/genotypes/support changed


class SampleChangeType(str, Enum):
    """How the carrier set of one variant entry changed."""

    CARRIER_ADDED = "CARRIER_ADDED"
    CARRIER_REMOVED = "CARRIER_REMOVED"
    GENOTYPE_CHANGED = "GENOTYPE_CHANGED"


@dataclass(frozen=True)
class SampleChange:
    """One sample-level difference on an entry between two versions."""

    change_type: SampleChangeType
    sample: str
    before: str | None = None   # genotype value ("./.", "0/0", "0/1", "1/1") before
    after: str | None = None    # genotype value after


@dataclass(frozen=True)
class EntryChange:
    """One entry-level difference between two callset versions."""

    change_type: EntryChangeType
    variant: str
    samples_before: tuple[str, ...] = ()
    samples_after: tuple[str, ...] = ()
    sample_changes: tuple[SampleChange, ...] = ()


@dataclass(frozen=True)
class ConflictChange:
    """A conflict that newly appeared or was resolved."""

    appeared: bool
    conflict_type: ConflictType
    variant: str
    sample: str | None
    detail: str


@dataclass(frozen=True)
class ChangeSummary:
    """Reproducible summary of differences between a baseline and a new version."""

    baseline_version: str | None
    version: str
    submitted_samples: tuple[str, ...]      # samples carried by this increment
    corrected_samples: tuple[str, ...]      # submitted samples that existed before
    new_samples: tuple[str, ...]            # submitted samples absent before
    entries_added: tuple[str, ...]
    entries_removed: tuple[str, ...]
    entry_changes: tuple[EntryChange, ...]
    conflicts_appeared: tuple[ConflictChange, ...]
    conflicts_resolved: tuple[ConflictChange, ...]


@dataclass(frozen=True)
class VersionedCallset:
    """One immutable callset version and the state it was derived from.

    Per-sample state (accepted normalized calls + sample-attributable
    rejections) lets later increments *replace* one sample wholesale instead
    of stacking new and old conclusions as extra evidence. ``batch_rejected``
    holds row-level structural rejections owned by one multi-sample batch
    (keyed by the frozen set of samples that batch declared).
    """

    version_id: str
    parent_version_id: str | None
    result: CallsetResult
    accepted_by_sample: dict[str, tuple[VariantRecord, ...]]
    rejected_by_sample: dict[str, tuple[Rejection, ...]]
    batch_rejected: dict[frozenset[str], tuple[Rejection, ...]]
    sample_order: tuple[str, ...]           # first-seen order, stable across versions
    delta_samples: frozenset[str] = frozenset()
    summary: ChangeSummary | None = None
