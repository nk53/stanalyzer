"""Strategy for independent or already chunked tasks."""

from typing import Callable, Iterable, TypeVar

from ..backend import ExecutionBackend


Task = TypeVar("Task")
Result = TypeVar("Result")


class IndependentStrategy:
    """Delegate independent tasks directly to the selected backend."""

    def execute(
        self,
        backend: ExecutionBackend,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        return backend.run(worker_fn, tasks, n_workers=n_workers)
