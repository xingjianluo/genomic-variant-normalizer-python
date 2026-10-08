"""HTTP API exposing the offline service.

Endpoints (all processing uses local data only):

* ``GET  /health``                  liveness / offline flag
* ``POST /normalize``               normalize one single-sample block
* ``POST /callset``                 normalize + merge a single/multi-sample block
* ``POST /increment``               apply a single/multi-sample increment onto a
                                    previously returned ``version_id``

Malformed input is never a transport error: responses are HTTP 200 with
``valid=false`` and a per-record ``rejected`` list carrying distinguishable
reason codes and exact locations. An empty submission or an empty sample name
on ``/normalize`` yields HTTP 422. An increment referencing an unknown/evicted
baseline yields 404 (``BASELINE_NOT_FOUND``); an increment carrying no samples
yields 422 (``EMPTY_INCREMENT``).
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from .models import CallsetResult, Conflict, Rejection
from .service import BaselineNotFound, CallsetService, EmptyIncrement

router = APIRouter()
service = CallsetService()


class SingleSampleRequest(BaseModel):
    text: str = Field(..., description="单样本变异文本，每行 CHROM POS REF ALT [GT]")
    sample: str = Field(..., min_length=1, description="样本名")


class MultiSampleRequest(BaseModel):
    text: str = Field(..., description="多样本变异文本，支持紧凑表头或 VCF 风格表头")


class IncrementRequest(BaseModel):
    text: str = Field(..., description="本次增量：仅含新增样本或需修正样本的文本")
    baseline_version: str = Field(..., min_length=1, description="上一轮返回的 version_id 基线标识")
    sample: str | None = Field(None, description="单样本形态时的样本名；多样本形态留空")


def _location_dict(loc) -> dict:
    return {"block": loc.block, "line_no": loc.line_no, "raw": loc.raw}


def serialize_rejection(r: Rejection) -> dict:
    return {
        "reason": r.reason.value,
        "detail": r.detail,
        "sample": r.sample,
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


def serialize_summary(summary) -> dict | None:
    if summary is None:
        return None

    def conflict_change(c) -> dict:
        return {
            "status": "appeared" if c.appeared else "resolved",
            "type": c.conflict_type.value,
            "variant": c.variant,
            "sample": c.sample,
            "detail": c.detail,
        }

    return {
        "baseline_version": summary.baseline_version,
        "version": summary.version,
        "submitted_samples": list(summary.submitted_samples),
        "new_samples": list(summary.new_samples),
        "corrected_samples": list(summary.corrected_samples),
        "entries_added": list(summary.entries_added),
        "entries_removed": list(summary.entries_removed),
        "entries_modified": [
            {
                "variant": c.variant,
                "samples_before": list(c.samples_before),
                "samples_after": list(c.samples_after),
                "sample_changes": [
                    {
                        "change": x.change_type.value,
                        "sample": x.sample,
                        "before": x.before,
                        "after": x.after,
                    }
                    for x in c.sample_changes
                ],
            }
            for c in summary.entry_changes
        ],
        "conflicts_appeared": [conflict_change(c) for c in summary.conflicts_appeared],
        "conflicts_resolved": [conflict_change(c) for c in summary.conflicts_resolved],
    }


def serialize_version(version) -> dict:
    body = serialize_result(version.result)
    body["version_id"] = version.version_id
    body["parent_version_id"] = version.parent_version_id
    body["change_summary"] = serialize_summary(version.summary)
    return body


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
        raise HTTPException(status_code=422, detail="text 与 sample 均不能为空")
    version = service.submit_single_version(req.text, sample)
    return serialize_version(version)


@router.post("/callset")
def build_callset(req: MultiSampleRequest):
    if not req.text.strip():
        raise HTTPException(status_code=422, detail="text 不能为空")
    # A multi-sample parser accepts single-sample text too (one inferred column),
    # unifying both input shapes under one endpoint.
    version = service.submit_multi_version(req.text)
    return serialize_version(version)


@router.post("/increment")
def increment_callset(req: IncrementRequest):
    if not req.text.strip():
        raise HTTPException(status_code=422, detail="text 不能为空")
    sample = req.sample.strip() if req.sample else None
    try:
        if sample:
            version = service.increment_single(req.text, sample, req.baseline_version)
        else:
            version = service.increment_multi(req.text, req.baseline_version)
    except BaselineNotFound:
        raise HTTPException(
            status_code=404,
            detail={
                "code": "BASELINE_NOT_FOUND",
                "baseline_version": req.baseline_version,
                "message": "基线版本不存在或已被淘汰，请基于最新 version_id 提交",
            },
        )
    except EmptyIncrement as exc:
        raise HTTPException(
            status_code=422,
            detail={"code": "EMPTY_INCREMENT", "message": str(exc)},
        )
    return serialize_version(version)
