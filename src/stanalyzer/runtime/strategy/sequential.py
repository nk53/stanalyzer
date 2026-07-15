"""Correctness-first strategy for ordered, stateful workloads."""

from typing import Callable, Iterable, TypeVar

from ..backend import ExecutionBackend, SequentialBackend


Task = TypeVar("Task")
Result = TypeVar("Result")


class SequentialStrategy:
    """Execute tasks locally in their supplied order."""

    def execute(
        self,
        backend: ExecutionBackend,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        del backend
        return SequentialBackend().run(
            worker_fn,
            tasks,
            n_workers=n_workers,
        )
