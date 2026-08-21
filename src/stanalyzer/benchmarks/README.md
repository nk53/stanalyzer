# ST-Analyzer benchmarks

`sta-benchmark` runs every configured analysis in a fresh process and records
wall time, CPU time, and peak RSS for the complete process tree (including
worker processes). It writes raw CSV, summary CSV, JSON metadata, and separate
stdout/stderr logs for each run.

Copy `example_manifest.json`, adjust the selections for your system, and run:

```console
sta-benchmark my_manifest.json --runs 3 --output-dir benchmark_results
```

The manifest uses argument arrays deliberately, so selections and paths do not
need shell quoting. `common_args` are placed before the analysis name; each
case's `args` are placed after it. Add one entry for every analysis applicable
to the dataset. A failed analysis is recorded and the rest continue unless
`--fail-fast` is supplied.

Peak RSS is sampled every 50 ms by default. For short analyses or sharper peak
detection, use `--sample-interval 0.01` (with slightly more sampling overhead).

Users can choose their input files without editing the manifest:

```console
sta-benchmark my_manifest.json --psf system.psf --traj run1.dcd run2.dcd
```

Use `--list-cases` to see the analyses configured by a manifest. Use one or
more `--analysis`/`-a` options to benchmark only selected cases:

```console
sta-benchmark my_manifest.json --list-cases
sta-benchmark my_manifest.json -a rmsd -a rmsf --psf system.psf --traj run.dcd
```

The selected files still need to be scientifically compatible with each
analysis. For example, a salt-bridge case needs charged residues and a
membrane-thickness case needs a membrane selection matching the system.
