#!/usr/bin/env bash
# Run test suite in the locked linux-64 pixi environment (after install_env.sh)
# pixi run activates env + sets tool data-dirs (LIBCIFPP_DATA_DIR etc.)
#
# Usage: run_test.sh [TEST_NAMES...]
#   No args = full suite; pass unittest paths to narrow, e.g.:
#     run_test.sh test_cli.CholTilt test_cli.CovAnalysis
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # <root>/src/stanalyzer/tests

export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0

# cd to tests/ for unittest imports (cwd on sys.path)
cd "$SCRIPT_DIR"

if [ "$#" -eq 0 ]; then
  pixi run test
else
  # test.sh doesn't pass args; invoke unittest directly
  pixi run python -m unittest -b "$@"
fi
echo TESTS_DONE