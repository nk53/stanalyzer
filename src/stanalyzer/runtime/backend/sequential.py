"""In-process execution backend."""

from typing import Callable, Iterable, TypeVar


Task = TypeVar("Task")
Result = TypeVar("Result")


class SequentialBackend:
    """Execute tasks serially in the calling process."""

    def run(
        self,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        del n_workers
        return [worker_fn(task) for task in tasks]

