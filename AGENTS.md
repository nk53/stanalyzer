# ST-Analyzer Agent Instructions

## Setup & dev commands

```bash
# Dev install (requires SETUPTOOLS_SCM_PRETEND_VERSION to avoid git tag failures)
./dev_install.sh                          # pip editable install
# or manually:
export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
pip install -e .

# Run CLI analysis
stanalyzer <analysis_name> [-h] [ARGS ...]  # run an analysis
stanalyzer -l                                # list all analyses
stanalyzer config                            # interactive project.json creator

# Run tests
./src/stanalyzer/tests/test.sh               # python -m unittest -b

# pixi workflow (official, recommended)
pixi install                                 # create/update project environment
pixi run test                                # run full test suite
pixi run smoke                               # CLI smoke check (stanalyzer -h && stanalyzer -l)

# Start dev web server (uvicorn, port 8000)
sta-server -r                                # -r enables hot-reload
```

Dev dependencies are in `dev_requirements.txt` (mypy, sphinx). Install with:
```bash
conda install -c conda-forge --file dev_requirements.txt
```

Linting and type checking:
- **Flake8**: per-file ignore on `src/stanalyzer/analysis/msd_membrane.py: E221` (config in `.flake8`)
- **Mypy**: `pyproject.toml [tool.mypy]` (check_untyped_defs, disallow_incomplete_defs, warn_* flags, import-untyped disabled; pydantic plugin)
- Tests use `unittest` — two files: `test_runtime.py`, `test_cli.py`

## Pull requests

Contributions must follow `CONTRIBUTING.md`. Read it before opening a pull
request. Every line of a contribution must be understood and explainable by
the human submitter — "the AI wrote it" is not an acceptable answer to a
review question. AI involvement must be documented per CONTRIBUTING.md.

## Architecture overview

Single package `stanalyzer`. Package layout:

| Directory | Purpose |
|---|---|
| `src/stanalyzer/cli/` | CLI entrypoints (`stanalyzer`, `sta-server`) and config handling |
| `src/stanalyzer/analysis/` | ~40 analysis modules. each must define `main(settings)` and `get_parser()` to be auto-discovered |
| `src/stanalyzer/workers/` | Background task workers (executor, contacts, pi_stacking, rdf, sasa, salt_bridge) |
| `src/stanalyzer/runtime/` | Scheduling/backends (sequential, process-based) and chunking |
| `src/stanalyzer/static/` | Web GUI assets (forms.js, style.css, jquery) |

Adding an analysis: create a `.py` file in `analysis/` with `main()` and `get_parser()`.

## Important constraints

- **Settings flow**: CLI reads `project.json` (default) for defaults. Web GUI writes this file. Arguments not set on CLI fall back to `project.json`; unset keys raise an error.
- **Trajectory files**: `stanalyzer` resolves globs sorted by **numeric** order (`\d+`). Braced glob syntax is supported.
- **Output paths**: Relative output paths are written under `<output_path>/<analysis_name>/`. Absolute paths bypass this.
- **Optional external deps**: secondary structure needs `dssp`, SASA needs `freesasa`, pore radius needs `hole2`. Install via conda-forge if needed.
- **Browser caching**: static files (`forms.js`, `style.css`) may be cached by the browser during dev. Use private mode or disable cache in dev tools.

## Build / packaging notes

- Version is managed by `setuptools_scm` (derived from git tags). For local development without tags, set `SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0`.
- Conda build files are in `conda.build/`.
- Entrypoints defined in `pyproject.toml`: `stanalyzer` → `stanalyzer.cli.stanalyzer:main`; `sta-server` → `stanalyzer.cli.stanalyzer:run_server`.
