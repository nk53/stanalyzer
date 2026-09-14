#!/usr/bin/env bash
# One-time setup: build the locked linux-64 pixi environment (inside docker)
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # <root>/src/stanalyzer/tests

cd "$REPO_ROOT"
# Cache pixi packages in workdir to avoid re-downloading in fresh containers
export PIXI_HOME="$REPO_ROOT/.pixi-home"
export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
pixi install --locked 2>&1 | tail -25

"$REPO_ROOT/.pixi/envs/default/bin/python" -c \
  "import numpy, scipy; print('numpy', numpy.__version__, '| scipy', scipy.__version__)"
echo INSTALL_DONE