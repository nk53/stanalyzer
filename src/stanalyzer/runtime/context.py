"""Facts about the environment in which an analysis will execute."""

from dataclasses import dataclass
import os
import sys
from typing import FrozenSet, Literal


Environment = Literal["desktop", "browser"]


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
    def detect_desktop(cls) -> "RuntimeContext":
        """Detect capabilities available to ordinary CPython."""
        return cls(
            environment="desktop",
            platform=sys.platform,
            logical_cpus=os.cpu_count() or 1,
            capabilities=frozenset({"process", "sequential"}),
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
