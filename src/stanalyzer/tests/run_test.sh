#!/usr/bin/env bash
# Run the test suite under the locked linux-64 pixi environment.
#
# Intended to run INSIDE the official pixi docker image
# (ghcr.io/prefix-dev/pixi:latest) after install_env.sh has built the
# environment (see README.md, "Dockerized linux-64 testing"). The script
# derives the repo root from its own location, so the mount point does not
# matter.
#
# `pixi run` activates the environment, so conda activation scripts run and
# tool data-dir env vars (e.g. LIBCIFPP_DATA_DIR for mkdssp) are set.
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

# cd to the tests dir so `python -m unittest` can import test_cli/test_runtime
# (cwd is on sys.path); pixi finds the project root by walking up from cwd.
cd "$SCRIPT_DIR"

if [ "$#" -eq 0 ]; then
  # Full suite: the pixi.toml `test` task (bash test.sh, cwd: tests dir).
  pixi run test
else
  # Narrow run: test.sh does not pass args through, so invoke unittest
  # directly via pixi run (which activates the environment).
  pixi run python -m unittest -b "$@"
fi
echo TESTS_DONE