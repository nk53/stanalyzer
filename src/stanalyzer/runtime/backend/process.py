"""Multiprocessing execution backend."""

from concurrent.futures import ProcessPoolExecutor
from typing import Callable, Iterable, TypeVar


Task = TypeVar("Task")
Result = TypeVar("Result")


class ProcessBackend:
    """Execute tasks in a local process pool."""

    def run(
        self,
        worker_fn: Callable[[Task], Result],
        tasks: Iterable[Task],
        *,
        n_workers: int,
    ) -> list[Result]:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            # Executor.map preserves input ordering and avoids exposing futures.
            return list(pool.map(worker_fn, tasks))

