import csv
import json
import sys

from stanalyzer.benchmarks.runner import (
    build_parser,
    load_cases,
    run_command,
    write_results,
)


def test_load_cases(tmp_path):
    manifest = tmp_path / "cases.json"
    manifest.write_text(json.dumps({"common_args": ["-s", "x.json"], "analyses": [{"name": "rmsd", "args": ["--sel", "all"]}]}))
    common, cases = load_cases(manifest)
    assert common == ["-s", "x.json"]
    assert cases[0]["name"] == "rmsd"


def test_run_command_and_results(tmp_path):
    row = run_command(
        "tiny",
        [sys.executable, "-c", "import time; x=bytearray(2_000_000); time.sleep(.1)"],
        1,
        tmp_path,
        0.01,
    )
    assert row["status"] == "success"
    assert row["wall_seconds"] >= 0
    assert row["peak_rss_mb"] > 0
    write_results([row], tmp_path)
    with (tmp_path / "benchmark_summary.csv").open() as stream:
        summary = list(csv.DictReader(stream))
    assert summary[0]["analysis"] == "tiny"
    assert summary[0]["successful_runs"] == "1"


def test_input_file_and_analysis_options():
    args = build_parser().parse_args([
        "cases.json", "--psf", "system.psf", "--traj", "a.dcd", "b.dcd",
        "-a", "rmsd", "-a", "rmsf",
    ])
    assert args.psf.name == "system.psf"
    assert [path.name for path in args.traj] == ["a.dcd", "b.dcd"]
    assert args.analysis == ["rmsd", "rmsf"]
