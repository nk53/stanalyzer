"""Public runtime executor facade."""

import os
from typing import Callable, Iterable, TypeVar

from .backend import ExecutionBackend, ProcessBackend, SequentialBackend
from .plan import ExecutionPlan
from .strategy import (
    ExecutionStrategy,
    IndependentStrategy,
    SequentialStrategy,
)


Task = TypeVar("Task")
Result = TypeVar("Result")


class RuntimeExecutor:
    """Execute analysis tasks through replaceable strategies and backends."""

    def __init__(
        self,
        n_workers: int | None = None,
        backend: ExecutionBackend | None = None,
        plan: ExecutionPlan | None = None,
        strategy: ExecutionStrategy | None = None,
    ) -> None:
        if plan is not None and n_workers is not None:
            raise ValueError("provide either plan or n_workers, not both")
        if plan is not None:
            n_workers = plan.n_workers
        if n_workers is None:
            n_workers = os.cpu_count() or 1
        if n_workers < 1:
            raise ValueError("n_workers must be at least 1")

        self.n_workers = n_workers
        if backend is not None:
            self.backend = backend
        elif plan is None or plan.backend == "process":
            self.backend = ProcessBackend()
        elif plan.backend == "sequential":
            self.backend = SequentialBackend()
        else:
            raise ValueError(f"unsupported backend in execution plan: {plan.backend}")

        if strategy is not None:
            self.strategy = strategy
        elif self.n_workers == 1:
            self.strategy = SequentialStrategy()
        elif plan is None or plan.strategy == "independent":
            self.strategy = IndependentStrategy()
        elif plan.strategy == "sequential":
            self.strategy = SequentialStrategy()
        else:
            raise ValueError(
                f"unsupported strategy in execution plan: {plan.strategy}"
            )

    def run(
        self,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
    ) -> list[Result]:
        """Return worker results in task order."""
        return self.strategy.execute(
            self.backend,
            worker_fn,
            tasks,
            n_workers=self.n_workers,
        )
