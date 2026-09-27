"""pytest 全局配置：让单测日志（输入摘要 / 判定依据）输出到控制台。"""
import logging
import os


def pytest_configure(config):
    level = logging.DEBUG if os.getenv("TEST_DEBUG") else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
