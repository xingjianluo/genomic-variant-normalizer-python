"""Shared pytest fixtures (skeleton)."""
import pytest


@pytest.fixture
def service():
    from variant_service.service import CallsetService
    return CallsetService()
