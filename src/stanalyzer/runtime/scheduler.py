"""Turn runtime facts and workload bounds into an execution plan."""

from .context import RuntimeContext
from .plan import ExecutionPlan


class RuntimeScheduler:
    """Create conservative plans without analysis-specific policy."""

    def __init__(self, context: RuntimeContext) -> None:
        self.context = context

    def create_plan(
        self,
        *,
        task_count: int | None = None,
        n_workers: int | None = None,
    ) -> ExecutionPlan:
        if task_count is not None and task_count < 1:
            raise ValueError("task_count must be at least 1")
        if n_workers is not None and n_workers < 1:
            raise ValueError("n_workers must be at least 1")

        workers = n_workers or self.context.logical_cpus
        workers = min(workers, self.context.logical_cpus)
        if task_count is not None:
            workers = min(workers, task_count)

        # A browser worker backend is deliberately not selected until it is
        # implemented and advertised as a context capability.
        backend = "process" if (
            workers > 1 and "process" in self.context.capabilities
        ) else "sequential"
        strategy = (
            "independent"
            if backend == "process"
            else "sequential"
        )
        return ExecutionPlan(
            backend=backend,
            n_workers=workers,
            strategy=strategy,
        )
