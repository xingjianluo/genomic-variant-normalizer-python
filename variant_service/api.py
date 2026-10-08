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

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from .models import CallsetResult, CallsetSnapshot, Conflict, Rejection
from .service import BaselineError, CallsetService

router = APIRouter()
service = CallsetService()


class SingleSampleRequest(BaseModel):
    text: str = Field(..., description="单样本变异文本，每行 CHROM POS REF ALT [GT]")
    sample: str = Field(..., min_length=1, description="样本名")


class MultiSampleRequest(BaseModel):
    text: str = Field(..., description="多样本变异文本，支持紧凑表头或 VCF 风格表头")


class IncrementSingleRequest(BaseModel):
    text: str = Field(..., description="新增/修正样本的单样本变异文本")
    sample: str = Field(..., min_length=1, description="样本名（该样本结论整体替换）")
    base_version: str = Field(..., min_length=1, description="上一轮返回的版本标识")


class IncrementMultiRequest(BaseModel):
    text: str = Field(..., description="新增/修正样本的多样本变异文本（两种表头均可）")
    base_version: str = Field(..., min_length=1, description="上一轮返回的版本标识")


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
        "version_id": result.version_id,
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
        raise HTTPException(status_code=422, detail="text 不能为空")
    # A multi-sample parser accepts single-sample text too (one inferred column),
    # unifying both input shapes under one endpoint.
    result = service.submit_multi(req.text)
    return serialize_result(result)


def serialize_entry_change(ch) -> dict:
    return {
        "variant": ch.variant.render(),
        "kind": ch.kind,
        "samples_added": list(ch.samples_added),
        "samples_removed": list(ch.samples_removed),
        "genotype_changes": [
            {"sample": s, "old": old, "new": new} for s, old, new in ch.genotype_changes
        ],
    }


def serialize_changes(snapshot: CallsetSnapshot) -> dict:
    cs = snapshot.changes
    return {
        "identical": cs.identical,
        "added_entries": [v.render() for v in cs.added_entries],
        "removed_entries": [v.render() for v in cs.removed_entries],
        "changed_entries": [serialize_entry_change(c) for c in cs.changed_entries],
        "samples_added": list(cs.samples_added),
        "samples_removed": list(cs.samples_removed),
        "new_conflicts": [serialize_conflict(c) for c in cs.new_conflicts],
        "resolved_conflicts": [serialize_conflict(c) for c in cs.resolved_conflicts],
        "new_rejections": [serialize_rejection(r) for r in cs.new_rejections],
        "resolved_rejections": [serialize_rejection(r) for r in cs.resolved_rejections],
    }


def serialize_snapshot(snapshot: CallsetSnapshot) -> dict:
    body = serialize_result(snapshot.result)
    body.update(
        {
            "version_id": snapshot.version_id,
            "parent_version": snapshot.parent_version,
            "root_version": snapshot.root_version,
            "sequence": snapshot.sequence,
            "samples": list(snapshot.samples),
            "touched_samples": list(snapshot.touched_samples),
            "changes": serialize_changes(snapshot),
        }
    )
    return body


def _raise_baseline(err: BaselineError) -> None:
    # Deterministic, distinguishable rejection: 422 for an unusable empty
    # increment, 404 for an unknown version, 409 for a stale/contended baseline.
    if err.error_type.value == "EMPTY_INCREMENT":
        raise HTTPException(
            status_code=422,
            detail={
                "reason": err.error_type.value,
                "detail": err.detail,
                "current_version": None,
                "conflicting_samples": [],
            },
        )
    status = 404 if err.error_type.value == "BASELINE_NOT_FOUND" else 409
    raise HTTPException(
        status_code=status,
        detail={
            "reason": err.error_type.value,
            "detail": err.detail,
            "current_version": err.current_version,
            "conflicting_samples": list(err.conflicting_samples),
        },
    )


@router.post("/increment/root/single")
def root_single(req: SingleSampleRequest):
    sample = req.sample.strip()
    if not sample or not req.text.strip():
        raise HTTPException(status_code=422, detail="text 与 sample 均不能为空")
    snapshot = service.submit_root_single(req.text, sample)
    return serialize_snapshot(snapshot)


@router.post("/increment/root/multi")
def root_multi(req: MultiSampleRequest):
    if not req.text.strip():
        raise HTTPException(status_code=422, detail="text 不能为空")
    snapshot = service.submit_root_multi(req.text)
    return serialize_snapshot(snapshot)


@router.post("/increment/single")
def increment_single(req: IncrementSingleRequest):
    sample = req.sample.strip()
    if not sample or not req.text.strip():
        raise HTTPException(status_code=422, detail="text 与 sample 均不能为空")
    try:
        snapshot = service.submit_increment_single(req.text, sample, req.base_version)
    except BaselineError as err:
        _raise_baseline(err)
    return serialize_snapshot(snapshot)


@router.post("/increment/multi")
def increment_multi(req: IncrementMultiRequest):
    if not req.text.strip():
        raise HTTPException(status_code=422, detail="text 不能为空")
    try:
        snapshot = service.submit_increment_multi(req.text, req.base_version)
    except BaselineError as err:
        _raise_baseline(err)
    return serialize_snapshot(snapshot)


@router.get("/version/{version_id}")
def get_version(version_id: str):
    snapshot = service.get_version(version_id)
    if snapshot is None:
        raise HTTPException(status_code=404, detail={"reason": "BASELINE_NOT_FOUND",
                                                     "detail": f"未知版本标识: {version_id}"})
    return serialize_snapshot(snapshot)
