"""FastAPI 路由：全程离线，仅使用本地解析/合并能力。"""
from __future__ import annotations

import logging
from typing import Any

from fastapi import APIRouter, HTTPException

from .models import Callset, SampleInput, SubmitResponse
from .service import VariantService

logger = logging.getLogger("variant_service.api")

router = APIRouter()
service = VariantService(max_concurrent=4)


def _summarize(callset: Callset) -> SubmitResponse:
    return SubmitResponse(
        callset_id=callset.callset_id,
        total_inputs=callset.total_inputs,
        accepted_calls=callset.accepted_calls,
        rejected=len(callset.rejected_records),
        variants=len(callset.variants),
    )


@router.get("/health")
def health() -> dict[str, Any]:
    return {"status": "ok", "offline": True, "cached_callsets": service.cached_count()}


@router.post("/callsets", response_model=SubmitResponse)
async def submit_callset(inputs: list[SampleInput]) -> SubmitResponse:
    """提交一批样本输入，返回 callset 概览（拒绝记录数而非丢弃）。"""
    if not inputs:
        raise HTTPException(status_code=400, detail="inputs 不能为空：至少提交一份样本输入")
    logger.info("收到提交: inputs=%d", len(inputs))
    callset = await service.submit(inputs)
    return _summarize(callset)


@router.post("/callsets/inspect", response_model=Callset)
async def inspect_callset(inputs: list[SampleInput]) -> Callset:
    """提交并返回完整 callset 明细（变异聚合、冲突标记、拒绝原因与位置）。"""
    if not inputs:
        raise HTTPException(status_code=400, detail="inputs 不能为空：至少提交一份样本输入")
    callset = await service.submit(inputs)
    return callset


@router.get("/callsets/{callset_id}", response_model=Callset)
def get_callset(callset_id: str) -> Callset:
    """按 ID 取回此前已生成的 callset（重复提交可直接 GET 复用）。"""
    callset = service.get_callset(callset_id)
    if callset is None:
        raise HTTPException(status_code=404, detail=f"未找到 callset: {callset_id}")
    return callset
