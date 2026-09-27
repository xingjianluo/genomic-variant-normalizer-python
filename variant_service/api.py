"""HTTP API exposing the offline service.

Endpoints (all processing uses local data only):

* ``GET  /health``                  liveness / offline flag
* ``POST /normalize``               normalize one single-sample block
* ``POST /callset``                 normalize + merge a single/multi-sample block

Malformed input is never a transport error: responses are HTTP 200 with
``valid=false`` and a per-record ``rejected`` list carrying distinguishable
reason codes and exact locations. An empty submission or an empty sample name
on ``/normalize`` yields HTTP 422.
"""
from __future__ import annotations

from fastapi import APIRouter
from pydantic import BaseModel, Field

from .models import CallsetResult, Conflict, Rejection
from .service import CallsetService

router = APIRouter()
service = CallsetService()


class SingleSampleRequest(BaseModel):
    text: str = Field(..., description="单样本变异文本，每行 CHROM POS REF ALT [GT]")
    sample: str = Field(..., min_length=1, description="样本名")


class MultiSampleRequest(BaseModel):
    text: str = Field(..., description="多样本变异文本，支持紧凑表头或 VCF 风格表头")


def _location_dict(loc) -> dict:
    return {"block": loc.block, "line_no": loc.line_no, "raw": loc.raw}


def serialize_rejection(r: Rejection) -> dict:
    return {
        "reason": r.reason.value,
        "detail": r.detail,
        "location": _location_dict(r.location),
    }


def serialize_conflict(c: Conflict) -> dict:
    return {
        "type": c.conflict_type.value,
        "variant": c.variant.render(),
        "sample": c.sample,
        "detail": c.detail,
        "locations": [_location_dict(l) for l in c.locations],
    }


def serialize_result(result: CallsetResult) -> dict:
    return {
        "content_hash": result.content_hash,
        "valid": len(result.rejected) == 0 and len(result.conflicts) == 0,
        "accepted_count": result.accepted_count,
        "entries": [
            {
                "variant": e.variant.render(),
                "chrom": e.variant.chrom,
                "pos": e.variant.pos,
                "ref": e.variant.ref,
                "alt": e.variant.alt,
                "samples": list(e.samples),
                "genotype_by_sample": {s: g.value for s, g in e.genotype_by_sample.items()},
                "support": [_location_dict(l) for l in e.support],
                "conflicts": [serialize_conflict(c) for c in e.conflicts],
            }
            for e in result.entries
        ],
        "rejected": [serialize_rejection(r) for r in result.rejected],
        "conflicts": [serialize_conflict(c) for c in result.conflicts],
        "summary": {
            "entry_count": len(result.entries),
            "rejected_count": len(result.rejected),
            "conflict_count": len(result.conflicts),
            "cache": service.cache_info(),
        },
    }


@router.get("/health")
def health():
    return {"status": "ok", "offline": True}


@router.post("/normalize")
def normalize_single(req: SingleSampleRequest):
    sample = req.sample.strip()
    if not sample or not req.text.strip():
        from fastapi import HTTPException
        raise HTTPException(status_code=422, detail="text 与 sample 均不能为空")
    result = service.submit_single(req.text, sample)
    return serialize_result(result)


@router.post("/callset")
def build_callset(req: MultiSampleRequest):
    if not req.text.strip():
        from fastapi import HTTPException
        raise HTTPException(status_code=422, detail="text 不能为空")
    # A multi-sample parser accepts single-sample text too (one inferred column),
    # unifying both input shapes under one endpoint.
    result = service.submit_multi(req.text)
    return serialize_result(result)
