"""Facts about the environment in which an analysis will execute."""

from dataclasses import dataclass
import os
import sys
from typing import FrozenSet, Literal


Environment = Literal["desktop", "browser"]

_BROWSER_PLATFORMS = frozenset({"emscripten", "wasi"})
_CPU_ALLOCATION_ENV_VARS = (
    "SLURM_CPUS_PER_TASK",
    "PBS_NP",
    "NSLOTS",
    "OMP_NUM_THREADS",
)


def _positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def _allocated_cpus() -> int:
    """Return the CPUs this process is actually allowed to use."""
    limits = [os.cpu_count() or 1]

    get_affinity = getattr(os, "sched_getaffinity", None)
    if get_affinity is not None:
        try:
            limits.append(len(get_affinity(0)))
        except (OSError, NotImplementedError):
            pass

    for variable in _CPU_ALLOCATION_ENV_VARS:
        value = _positive_int(os.environ.get(variable))
        if value is not None:
            limits.append(value)

    return max(1, min(limits))


def _available_memory_bytes() -> int | None:
    """Best-effort available-memory detection without extra dependencies."""
    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, OSError, ValueError):
        return None

    if not isinstance(pages, int) or not isinstance(page_size, int):
        return None
    if pages < 1 or page_size < 1:
        return None
    return pages * page_size


@dataclass(frozen=True)
class RuntimeContext:
    """Platform capabilities supplied to the scheduler.

    Browser values are injected by the JavaScript host. Keeping detection out
    of the scheduler lets the same planning logic run on desktop and Pyodide.
    """

    environment: Environment
    platform: str
    logical_cpus: int
    capabilities: FrozenSet[str] = frozenset()
    available_memory_bytes: int | None = None

    def __post_init__(self) -> None:
        if self.logical_cpus < 1:
            raise ValueError("logical_cpus must be at least 1")
        if self.available_memory_bytes is not None \
                and self.available_memory_bytes < 1:
            raise ValueError("available_memory_bytes must be at least 1")

    @classmethod
    def detect(cls) -> "RuntimeContext":
        """Detect a safe runtime context for CPython or WebAssembly."""
        if sys.platform in _BROWSER_PLATFORMS:
            return cls.from_browser(
                logical_cpus=os.cpu_count() or 1,
            )
        return cls.detect_desktop()

    @classmethod
    def detect_desktop(cls) -> "RuntimeContext":
        """Detect capabilities and resource limits for ordinary CPython."""
        return cls(
            environment="desktop",
            platform=sys.platform,
            logical_cpus=_allocated_cpus(),
            capabilities=frozenset({"process", "sequential"}),
            available_memory_bytes=_available_memory_bytes(),
        )

    @classmethod
    def from_browser(
        cls,
        *,
        logical_cpus: int,
        shared_array_buffer: bool = False,
        available_memory_bytes: int | None = None,
    ) -> "RuntimeContext":
        """Create a context from browser values such as hardwareConcurrency."""
        capabilities = {"sequential", "webassembly"}
        if shared_array_buffer:
            capabilities.add("shared-array-buffer")
        return cls(
            environment="browser",
            platform="browser",
            logical_cpus=logical_cpus,
            capabilities=frozenset(capabilities),
            available_memory_bytes=available_memory_bytes,
        )
