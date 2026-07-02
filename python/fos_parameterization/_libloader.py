"""Locate and load libfos_parameterization.so for ctypes."""
from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import Iterator, Optional

_LIB_NAME = "libfos_parameterization.so"


def _candidate_paths() -> Iterator[Path]:
    env = os.environ.get("FOS_PARAM_LIB")
    if env:
        yield Path(env)
    # Dev fallback: this file is <repo>/python/fos_parameterization/_libloader.py
    repo_root = Path(__file__).resolve().parents[2]
    yield repo_root / "build" / "release" / _LIB_NAME
    yield repo_root / "cmake-build-release" / _LIB_NAME
    yield repo_root / "cmake-build-debug" / _LIB_NAME
    # System search path
    yield Path(_LIB_NAME)


def load_library() -> ctypes.CDLL:
    """Load the shared library, trying FOS_PARAM_LIB, build dirs, then the system path."""
    last_error: Optional[OSError] = None
    for path in _candidate_paths():
        try:
            return ctypes.CDLL(str(path))
        except OSError as exc:
            last_error = exc
    raise OSError(
        f"Could not load {_LIB_NAME}. Build it (cmake --build build/release) or set "
        f"FOS_PARAM_LIB to its full path. Last error: {last_error}"
    )
