# ST-Analyzer benchmark runner

Runs configured analyses in separate processes and reports:

- **Runtime:** elapsed wall-clock time.
- **Peak RSS:** sampled memory for the analysis and its workers.

Run from a directory containing `project.json`. The selected Python environment
must contain ST-Analyzer and its analysis dependencies.

```console
sta-benchmark example_manifest.json \
  --psf /path/to/system.psf \
  --traj /path/to/run1.dcd /path/to/run2.dcd \
  --runs 3 \
  --output-dir benchmark_results \
  --python /path/to/environment/bin/python
```

List cases or run selected analyses:

```console
sta-benchmark example_manifest.json --list-cases
sta-benchmark example_manifest.json -a rmsd -a rmsf --runs 3
```

Results are written to:

```text
benchmark_results/benchmark_summary.csv
benchmark_results/benchmark_raw.csv
benchmark_results/benchmark_results.json
benchmark_results/logs/
```

The example manifest is a starter set; its analyses and selections must match
the supplied dataset. Failed cases are recorded while remaining cases continue.

Memory is sampled every 50 ms. Brief peaks may be missed, and shared worker
memory may be counted more than once. Compare results from the same machine and
environment.
