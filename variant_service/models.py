"""统一记录模型与枚举（骨架：仅类型定义）。"""
from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field


class InputFormat(str, Enum):
    SINGLE = "single"  # sample_id<TAB>chrom<TAB>pos<TAB>ref<TAB>alt[<TAB>gt]
    VCF = "vcf"        # 多样本 VCF 文本


class Genotype(str, Enum):
    HOM_REF = "HOM_REF"
    HET = "HET"
    HOM_ALT = "HOM_ALT"
    UNKNOWN = "UNKNOWN"


class RejectReason(str, Enum):
    MALFORMED_LINE = "MALFORMED_LINE"
    WRONG_COLUMN_COUNT = "WRONG_COLUMN_COUNT"
    SAMPLE_COLUMN_MISMATCH = "SAMPLE_COLUMN_MISMATCH"
    MISSING_FIELD = "MISSING_FIELD"
    INVALID_POSITION = "INVALID_POSITION"
    INVALID_ALLELE = "INVALID_ALLELE"
    INVALID_CHROM = "INVALID_CHROM"
    INVALID_GENOTYPE = "INVALID_GENOTYPE"
    DUPLICATE_SAMPLE_HEADER = "DUPLICATE_SAMPLE_HEADER"
    EMPTY_INPUT = "EMPTY_INPUT"


class VariantCall(BaseModel):
    sample_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: Genotype
    source_line: int


class RejectedRecord(BaseModel):
    source: str = ""
    line: int
    reason: RejectReason
    detail: str
    raw: str = ""


class SampleInput(BaseModel):
    sample_id: str = ""
    format: InputFormat
    text: str


class SubmitResponse(BaseModel):
    callset_id: str
    total_inputs: int
    accepted_calls: int
    rejected: int
    variants: int


class CallsetVariant(BaseModel):
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype_summary: dict[str, list[str]] = Field(default_factory=dict)
    samples: list[str] = Field(default_factory=list)
    call_count: int = 0
    conflicting: bool = False
    conflict_types: list[str] = Field(default_factory=list)
    conflicts: list[dict] = Field(default_factory=list)


class Callset(BaseModel):
    callset_id: str
    total_inputs: int
    accepted_calls: int
    rejected_records: list[RejectedRecord] = Field(default_factory=list)
    variants: list[CallsetVariant] = Field(default_factory=list)


class ParseResult(BaseModel):
    calls: list[VariantCall] = Field(default_factory=list)
    rejected: list[RejectedRecord] = Field(default_factory=list)
