"""Contracts implemented by execution backends."""

from typing import Callable, Iterable, Protocol, TypeVar


Task = TypeVar("Task")
Result = TypeVar("Result")


class ExecutionBackend(Protocol):
    """Run independent tasks without exposing the execution mechanism."""

    def run(
        self,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        """Return results in the same order as the supplied tasks."""

