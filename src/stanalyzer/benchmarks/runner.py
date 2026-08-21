"""Run ST-Analyzer analyses in isolation and measure time and peak memory.

Each benchmark case is a normal ``stanalyzer`` command.  Running cases in
separate processes makes the peak-memory value independent for every run.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import resource
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


FIELDS = (
    "analysis", "run", "status", "returncode", "wall_seconds",
    "cpu_seconds", "peak_rss_mb", "command", "stdout", "stderr",
)

CLI_CODE = "from stanalyzer.cli.stanalyzer import main; main()"


def _rss_process_tree_kb(root_pid: int) -> int:
    """Return summed RSS for root_pid and its descendants (best effort)."""
    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="], capture_output=True,
            text=True, check=False, timeout=5,
        )
        rows = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    except (OSError, subprocess.SubprocessError, ValueError):
        return 0

    children: dict[int, list[int]] = {}
    rss: dict[int, int] = {}
    for pid, ppid, rss_kb in rows:
        children.setdefault(ppid, []).append(pid)
        rss[pid] = rss_kb
    pending = [root_pid]
    descendants: set[int] = set()
    while pending:
        pid = pending.pop()
        if pid in descendants:
            continue
        descendants.add(pid)
        pending.extend(children.get(pid, ()))
    return sum(rss.get(pid, 0) for pid in descendants)


def _resource_rss_kb(value: float) -> float:
    """Normalize ru_maxrss (bytes on macOS, KiB on Linux) to KiB."""
    return value / 1024 if sys.platform == "darwin" else value


def run_command(
    analysis: str,
    command: list[str],
    run_id: int,
    output_dir: Path,
    sample_interval: float = 0.05,
) -> dict[str, Any]:
    """Benchmark one command, including any worker child processes."""
    case_dir = output_dir / "logs" / analysis
    case_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = case_dir / f"run-{run_id}.stdout.log"
    stderr_path = case_dir / f"run-{run_id}.stderr.log"
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    peak_kb = 0

    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(
            command, stdout=stdout, stderr=stderr, stdin=subprocess.DEVNULL,
            start_new_session=True,
        )
        while True:
            peak_kb = max(peak_kb, _rss_process_tree_kb(process.pid))
            if process.poll() is not None:
                break
            time.sleep(sample_interval)
        returncode = process.returncode

    wall = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime + after.ru_stime) - (before.ru_utime + before.ru_stime)
    # Some restricted environments disallow `ps`. ru_maxrss is a useful
    # fallback for the direct analysis process, though it cannot sum workers.
    if peak_kb == 0:
        peak_kb = int(_resource_rss_kb(after.ru_maxrss))
    return {
        "analysis": analysis,
        "run": run_id,
        "status": "success" if returncode == 0 else "failed",
        "returncode": returncode,
        "wall_seconds": round(wall, 6),
        "cpu_seconds": round(cpu, 6),
        "peak_rss_mb": round(peak_kb / 1024, 3),
        "command": shlex.join(command),
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def load_cases(path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    """Load common CLI arguments and analysis cases from a JSON manifest."""
    with path.open() as stream:
        data = json.load(stream)
    common = data.get("common_args", [])
    cases = data.get("analyses")
    if not isinstance(common, list) or not all(isinstance(x, str) for x in common):
        raise ValueError("common_args must be a list of strings")
    if not isinstance(cases, list) or not cases:
        raise ValueError("analyses must be a non-empty list")
    for case in cases:
        if not isinstance(case, dict) or not isinstance(case.get("name"), str):
            raise ValueError("each analysis needs a string name")
        if not isinstance(case.get("args", []), list):
            raise ValueError(f"args for {case.get('name', '?')} must be a list")
    return common, cases


def write_results(
    rows: list[dict[str, Any]], output_dir: Path,
    sample_interval: float | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "benchmark_raw.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for name in sorted({row["analysis"] for row in rows}):
        group = [row for row in rows if row["analysis"] == name]
        ok = [row for row in group if row["status"] == "success"]
        summary.append({
            "analysis": name,
            "successful_runs": len(ok),
            "failed_runs": len(group) - len(ok),
            "wall_mean_seconds": round(statistics.mean(r["wall_seconds"] for r in ok), 6) if ok else "",
            "wall_min_seconds": round(min(r["wall_seconds"] for r in ok), 6) if ok else "",
            "peak_rss_max_mb": round(max(r["peak_rss_mb"] for r in ok), 3) if ok else "",
            "cpu_mean_seconds": round(statistics.mean(r["cpu_seconds"] for r in ok), 6) if ok else "",
        })
    fields = tuple(summary[0]) if summary else ()
    with (output_dir / "benchmark_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)
    metadata = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "platform": platform.platform(),
        "python": sys.version,
        "sample_interval_seconds": sample_interval,
        "results": rows,
    }
    (output_dir / "benchmark_results.json").write_text(json.dumps(metadata, indent=2) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark ST-Analyzer runtime and peak process-tree memory")
    parser.add_argument("manifest", type=Path, help="JSON benchmark manifest")
    parser.add_argument("-n", "--runs", type=int, default=3, help="runs per analysis (default: 3)")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("benchmark_results"))
    parser.add_argument("--sample-interval", type=float, default=0.05, metavar="SECONDS")
    parser.add_argument("--python", default=sys.executable, help="Python executable used to launch stanalyzer")
    parser.add_argument("--psf", type=Path, help="override the topology file for every selected analysis")
    parser.add_argument(
        "--traj", type=Path, nargs="+", metavar="FILE",
        help="override trajectory file(s) for every selected analysis",
    )
    parser.add_argument(
        "-a", "--analysis", action="append", metavar="NAME",
        help="run only this manifest case; repeat for multiple analyses",
    )
    parser.add_argument(
        "--list-cases", action="store_true",
        help="list analyses configured in the manifest and exit",
    )
    parser.add_argument("--fail-fast", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.runs < 1 or args.sample_interval <= 0:
        raise SystemExit("--runs and --sample-interval must be positive")
    common, cases = load_cases(args.manifest)
    if args.list_cases:
        print("\n".join(case["name"] for case in cases))
        return 0
    if args.analysis:
        requested = set(args.analysis)
        known = {case["name"] for case in cases}
        unknown = sorted(requested - known)
        if unknown:
            raise SystemExit(
                "analysis not configured in manifest: " + ", ".join(unknown)
                + "; use --list-cases"
            )
        cases = [case for case in cases if case["name"] in requested]

    input_args: list[str] = []
    if args.psf:
        topology = args.psf.expanduser().resolve()
        if not topology.is_file():
            raise SystemExit(f"topology file does not exist: {topology}")
        input_args.extend(["--psf", str(topology)])
    if args.traj:
        trajectories = [path.expanduser().resolve() for path in args.traj]
        missing = [str(path) for path in trajectories if not path.is_file()]
        if missing:
            raise SystemExit("trajectory file does not exist: " + ", ".join(missing))
        input_args.append("--traj")
        input_args.extend(str(path) for path in trajectories)

    rows: list[dict[str, Any]] = []
    for case in cases:
        name = case["name"]
        for run_id in range(1, args.runs + 1):
            # Do not use ``python -m stanalyzer.cli.stanalyzer`` here. The
            # package's cli/__init__.py imports that module, so runpy would
            # execute a second copy as __main__ and split its module globals.
            command = [
                args.python, "-c", CLI_CODE, *common, name,
                *input_args, *case.get("args", []),
            ]
            print(f"[{name}] run {run_id}/{args.runs}", flush=True)
            row = run_command(name, command, run_id, args.output_dir, args.sample_interval)
            rows.append(row)
            print(f"  {row['status']}: {row['wall_seconds']:.3f}s, peak RSS {row['peak_rss_mb']:.1f} MB", flush=True)
            if row["status"] == "failed" and args.fail_fast:
                write_results(rows, args.output_dir, args.sample_interval)
                return 1
    write_results(rows, args.output_dir, args.sample_interval)
    print(f"Results: {args.output_dir / 'benchmark_summary.csv'}")
    return 1 if any(row["status"] == "failed" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
