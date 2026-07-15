"""Contract implemented by runtime execution strategies."""

from typing import Callable, Iterable, Protocol, TypeVar

from ..backend import ExecutionBackend


Task = TypeVar("Task")
Result = TypeVar("Result")


class ExecutionStrategy(Protocol):
    """Coordinate task execution without entering worker hot loops."""

    def execute(
        self,
        backend: ExecutionBackend,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        """Execute tasks according to the strategy's dependency behavior."""
