"""service 测试：确定性 ID、重复提交复用、并发不污染、万行性能。"""
import asyncio
import logging
import time

from variant_service.models import SampleInput
from variant_service.service import VariantService, compute_callset_id

logger = logging.getLogger("test.service")


def _single_input(sid, lines):
    return SampleInput(sample_id=sid, format="single", text="\n".join(lines))


def test_deterministic_id_independent_of_order():
    a = [
        _single_input("S1", ["chr1\t10\tA\tG"]),
        _single_input("S2", ["chr1\t20\tC\tT"]),
    ]
    b = list(reversed(a))
    id_a, id_b = compute_callset_id(a), compute_callset_id(b)
    logger.info("输入摘要 两份顺序互换的输入 -> id=%s / %s 判定一致=%s", id_a, id_b, id_a == id_b)
    assert id_a == id_b


async def _concurrent():
    svc = VariantService(max_concurrent=8)
    inp = [_single_input("S1", ["chr1\t10\tA\tG\tHET", "chr1\t20\tC\tT\tHOM_ALT"])]

    # 10 个相同请求并发
    results = await asyncio.gather(*[svc.submit(inp) for _ in range(10)])
    ids = {r.callset_id for r in results}
    assert len(ids) == 1
    assert svc.cached_count() == 1
    logger.info("输入摘要 10 个并发重复请求 -> 唯一 callset=%s, 缓存条目=%d", ids.pop(), svc.cached_count())

    # 重叠样本集合不互相污染
    shared = "chr1\t10\tA\tG\tHET"
    inp_a = [_single_input("S1", [shared]), _single_input("S2", ["chr1\t20\tC\tT\tHET"])]
    inp_b = [_single_input("S1", [shared]), _single_input("S3", ["chr2\t30\tG\tA\tHET"])]
    ra, rb = await asyncio.gather(svc.submit(inp_a), svc.submit(inp_b))
    assert ra.callset_id != rb.callset_id
    samples_a = {s for v in ra.variants for s in v.samples}
    samples_b = {s for v in rb.variants for s in v.samples}
    logger.info("判定依据: 重叠输入 A 样本=%s B 样本=%s，无交叉污染=%s",
                samples_a, samples_b, samples_a == {"S1", "S2"} and samples_b == {"S1", "S3"})
    assert samples_a == {"S1", "S2"}
    assert samples_b == {"S1", "S3"}
    # 重复提交返回同一对象（复用）
    ra2 = await svc.submit(inp_a)
    assert ra2 is ra
    return svc


def test_concurrent_and_overlapping_inputs():
    asyncio.run(_concurrent())


def test_repeated_submission_stable():
    svc = VariantService()
    inp = [_single_input("S1", ["chr1\t10\tA\tG\tHET", "chr1\t10\tA\tG\tHET"])]
    r1 = asyncio.run(svc.submit(inp))
    r2 = asyncio.run(svc.submit(inp))
    assert r1.model_dump() == r2.model_dump()
    assert r1.accepted_calls == 1  # 重复不重复计数
    logger.info("判定: 重复行去重 accepted=%d，两次提交结论完全一致", r1.accepted_calls)


def test_ten_thousand_lines_performance():
    lines = [f"chr1\t{100 + i}\tA\tG\tHET" for i in range(10_000)]
    inp = [_single_input("S1", lines)]
    svc = VariantService(max_concurrent=4)

    async def go():
        # 重请求 + 一个轻请求并发，轻请求不应被长时间饿死
        heavy = asyncio.create_task(svc.submit(inp))
        await asyncio.sleep(0)
        t0 = time.perf_counter()
        light = await svc.submit([_single_input("S9", ["chr9\t1\tC\tT"])])
        light_ms = (time.perf_counter() - t0) * 1000
        heavy_result = await heavy
        return heavy_result, light, light_ms

    t0 = time.perf_counter()
    heavy_result, light, light_ms = asyncio.run(go())
    elapsed = time.perf_counter() - t0
    logger.info("输入摘要 10000 行 -> variants=%d 耗时=%.2fs；并发轻请求耗时=%.1fms",
                len(heavy_result.variants), elapsed, light_ms)
    assert len(heavy_result.variants) == 10_000
    assert len(light.variants) == 1
    assert elapsed < 15.0
