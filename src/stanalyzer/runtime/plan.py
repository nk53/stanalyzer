"""Immutable decisions produced by the runtime scheduler."""

from dataclasses import dataclass


@dataclass(frozen=True)
class ExecutionPlan:
    """Execution policy consumed by RuntimeExecutor."""

    backend: str
    n_workers: int
    strategy: str = "independent"
    chunk_size: int | None = None

    def __post_init__(self) -> None:
        if not self.backend:
            raise ValueError("backend must not be empty")
        if not self.strategy:
            raise ValueError("strategy must not be empty")
        if self.n_workers < 1:
            raise ValueError("n_workers must be at least 1")
        if self.chunk_size is not None and self.chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")
