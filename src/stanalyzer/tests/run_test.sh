#!/usr/bin/env bash
# Run the test suite under the locked linux-64 pixi environment.
#
# Intended to run INSIDE a docker container after install_env.sh has built
# the environment (see README.md, "Dockerized linux-64 testing"). The script
# derives the repo root from its own location, so the mount point does not
# matter.
#
# Usage: run_test.sh [TEST_NAMES...]
#   With no arguments, runs the full suite - identical to `pixi run test`
#   (which is `bash test.sh` = `python -m unittest -b`, cwd: tests dir).
#   Pass unittest paths to narrow, e.g.:
#     run_test.sh test_cli.CholTilt test_cli.CovAnalysis
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # <root>/src/stanalyzer/tests

export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
export PATH="$REPO_ROOT/.pixi/envs/default/bin:$PATH"

cd "$SCRIPT_DIR"
python -m unittest -b "$@"
echo TESTS_DONE