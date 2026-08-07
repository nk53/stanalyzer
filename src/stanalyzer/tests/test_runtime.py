import unittest
from unittest.mock import patch

from stanalyzer.workers.executor import FrameExecutor
from stanalyzer.runtime import (
    ExecutionPlan,
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    SequentialStrategy,
    chunk_frames,
)


def square(value: int) -> int:
    return value * value


class RecordingBackend:
    def __init__(self) -> None:
        self.n_workers = None

    def run(self, worker_fn, tasks, *, n_workers):
        self.n_workers = n_workers
        return [worker_fn(task) for task in tasks]


class RuntimeExecutorTest(unittest.TestCase):
    def test_rejects_invalid_worker_count(self):
        with self.assertRaisesRegex(ValueError, "at least 1"):
            RuntimeExecutor(n_workers=0)

    def test_one_worker_executes_sequentially(self):
        backend = RecordingBackend()
        executor = RuntimeExecutor(n_workers=1, backend=backend)

        self.assertEqual(executor.run(square, [3, 1, 2]), [9, 1, 4])
        self.assertIsNone(backend.n_workers)

    def test_delegates_to_configured_backend(self):
        backend = RecordingBackend()
        executor = RuntimeExecutor(n_workers=3, backend=backend)

        self.assertEqual(executor.run(square, [3, 1, 2]), [9, 1, 4])
        self.assertEqual(backend.n_workers, 3)

    def test_process_backend_preserves_task_order(self):
        executor = RuntimeExecutor(n_workers=2)

        self.assertEqual(executor.run(square, [3, 1, 2]), [9, 1, 4])

    def test_frame_executor_remains_compatible(self):
        executor = FrameExecutor(n_workers=1)

        self.assertEqual(executor.run(square, [2, 4]), [4, 16])

    def test_execution_plan_configures_executor(self):
        plan = ExecutionPlan(backend="sequential", n_workers=3)

        self.assertEqual(
            RuntimeExecutor(plan=plan).run(square, [2, 4]),
            [4, 16],
        )

    def test_sequential_strategy_bypasses_parallel_backend(self):
        backend = RecordingBackend()
        executor = RuntimeExecutor(
            n_workers=4,
            backend=backend,
            strategy=SequentialStrategy(),
        )

        self.assertEqual(executor.run(square, [3, 1, 2]), [9, 1, 4])
        self.assertIsNone(backend.n_workers)

    def test_rejects_unknown_plan_strategy(self):
        plan = ExecutionPlan(
            backend="process",
            n_workers=2,
            strategy="unknown",
        )

        with self.assertRaisesRegex(ValueError, "unsupported strategy"):
            RuntimeExecutor(plan=plan)


class RuntimeSchedulerTest(unittest.TestCase):
    def test_desktop_plan_is_bounded_by_tasks_and_cpus(self):
        context = RuntimeContext(
            environment="desktop",
            platform="test",
            logical_cpus=8,
            capabilities=frozenset({"process", "sequential"}),
        )

        plan = RuntimeScheduler(context).create_plan(task_count=3)

        self.assertEqual(plan.n_workers, 3)
        self.assertEqual(plan.backend, "process")
        self.assertEqual(plan.strategy, "independent")

    def test_browser_context_does_not_claim_unimplemented_backend(self):
        context = RuntimeContext.from_browser(
            logical_cpus=12,
            shared_array_buffer=True,
        )

        plan = RuntimeScheduler(context).create_plan(task_count=20)

        self.assertEqual(plan.n_workers, 1)
        self.assertEqual(plan.backend, "sequential")
        self.assertEqual(plan.strategy, "sequential")
        self.assertIn("webassembly", context.capabilities)
        self.assertIn("shared-array-buffer", context.capabilities)

    def test_auto_worker_count_is_conservative(self):
        context = RuntimeContext(
            environment="desktop",
            platform="test",
            logical_cpus=32,
            capabilities=frozenset({"process", "sequential"}),
            available_memory_bytes=64 * 1024 ** 3,
        )

        plan = RuntimeScheduler(context).create_plan(task_count=100)

        self.assertEqual(plan.n_workers, 4)
        self.assertEqual(plan.backend, "process")

    def test_auto_worker_count_respects_available_memory(self):
        context = RuntimeContext(
            environment="desktop",
            platform="test",
            logical_cpus=16,
            capabilities=frozenset({"process", "sequential"}),
            available_memory_bytes=2 * 1024 ** 3,
        )

        plan = RuntimeScheduler(context).create_plan(task_count=100)

        self.assertEqual(plan.n_workers, 2)

    def test_explicit_worker_count_is_only_resource_bounded(self):
        context = RuntimeContext(
            environment="desktop",
            platform="test",
            logical_cpus=8,
            capabilities=frozenset({"process", "sequential"}),
            available_memory_bytes=1024 ** 3,
        )

        plan = RuntimeScheduler(context).create_plan(
            task_count=100,
            n_workers=6,
        )

        self.assertEqual(plan.n_workers, 6)

    def test_detect_uses_browser_context_on_webassembly(self):
        with patch("stanalyzer.runtime.context.sys.platform", "emscripten"):
            with patch("stanalyzer.runtime.context.os.cpu_count", return_value=8):
                context = RuntimeContext.detect()

        self.assertEqual(context.environment, "browser")
        self.assertNotIn("process", context.capabilities)


class ChunkFramesTest(unittest.TestCase):
    def test_splits_frames_evenly_without_losing_frames(self):
        self.assertEqual(
            chunk_frames(n_frames=10, n_workers=3),
            [(0, 4), (4, 8), (8, 10)],
        )

    def test_empty_trajectory_has_no_chunks(self):
        self.assertEqual(chunk_frames(n_frames=0, n_workers=4), [])

    def test_rejects_invalid_values(self):
        with self.assertRaisesRegex(ValueError, "non-negative"):
            chunk_frames(n_frames=-1, n_workers=1)
        with self.assertRaisesRegex(ValueError, "at least 1"):
            chunk_frames(n_frames=1, n_workers=0)


if __name__ == "__main__":
    unittest.main()
