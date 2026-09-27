"""离线处理服务：解析、合并、缓存与并发控制。

设计要点：
* callset_id 由归一化后的输入集合（sample_id/format/text，排序后）做 SHA-256
  得到，因此同一批输入无论请求顺序、无论提交多少次，ID 与结论都一致。
* 每个 callset_id 一把锁：并发提交相同/重叠输入时只计算一次，其余请求复用；
  不同输入并行计算，互不阻塞、互不污染（全程无共享可变工作状态）。
* 解析与合并是 CPU 密集操作，放到有界 ThreadPoolExecutor，避免阻塞事件循环；
  线程池大小即为同时进行的重计算上限，防止万行级批量请求拖垮其它请求。
"""
from __future__ import annotations

import asyncio
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor

from .merge import merge_calls
from .models import Callset, SampleInput
from .normalization import parse_input


def compute_callset_id(inputs: list[SampleInput]) -> str:
    payload = sorted(
        (obj.sample_id, obj.format.value if hasattr(obj.format, "value") else str(obj.format), obj.text)
        for obj in inputs
    )
    blob = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def _process(inputs: list[SampleInput]) -> Callset:
    """纯函数：解析全部输入 -> 合并。线程池中执行，无共享状态。"""
    callset_id = compute_callset_id(inputs)
    all_calls = []
    all_rejected = []
    for obj in inputs:
        fmt = obj.format.value if hasattr(obj.format, "value") else str(obj.format)
        result = parse_input(obj.sample_id, fmt, obj.text)
        all_calls.extend(result.calls)
        all_rejected.extend(result.rejected)
    return merge_calls(callset_id, len(inputs), all_calls, all_rejected)


class VariantService:
    def __init__(self, max_concurrent: int = 4) -> None:
        self._cache: dict[str, Callset] = {}
        self._key_locks: dict[str, asyncio.Lock] = {}
        self._guard = asyncio.Lock()
        self._executor = ThreadPoolExecutor(
            max_workers=max_concurrent, thread_name_prefix="variant-worker"
        )

    async def submit(self, inputs: list[SampleInput]) -> Callset:
        callset_id = compute_callset_id(inputs)

        # 快速路径：已有结论直接复用
        cached = self._cache.get(callset_id)
        if cached is not None:
            return cached

        async with self._guard:
            lock = self._key_locks.setdefault(callset_id, asyncio.Lock())

        async with lock:
            cached = self._cache.get(callset_id)
            if cached is not None:
                return cached  # 并发重复提交：复用前一个请求的结论
            loop = asyncio.get_running_loop()
            callset = await loop.run_in_executor(self._executor, _process, inputs)
            self._cache[callset_id] = callset
            return callset

    def get_callset(self, callset_id: str) -> Callset | None:
        return self._cache.get(callset_id)

    def cached_count(self) -> int:
        return len(self._cache)
