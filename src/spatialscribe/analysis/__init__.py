"""Analysis engine: pure functions shared by the wizard UI and the Claude copilot."""

from .io import SpatialSample, load, detect_platform, build_control_mask
from .backend import get_backend, gpu_available

__all__ = [
    "SpatialSample",
    "load",
    "detect_platform",
    "build_control_mask",
    "get_backend",
    "gpu_available",
]
