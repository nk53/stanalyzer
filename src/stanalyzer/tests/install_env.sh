#!/usr/bin/env bash
# One-time setup: build the locked linux-64 environment.
#
# Intended to run INSIDE the official pixi docker image
# (ghcr.io/prefix-dev/pixi:latest), which already provides the pixi binary,
# with a copy of the repo mounted somewhere (see README.md, "Dockerized
# linux-64 testing"). The script derives the repo root from its own location,
# so the mount point does not matter.
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # <root>/src/stanalyzer/tests

cd "$REPO_ROOT"
# Persist the pixi package cache inside the mounted workdir so the env does
# not need to be re-downloaded the next time install_env.sh is run in a
# fresh container (the default cache lives in the container's $HOME).
export PIXI_HOME="$REPO_ROOT/.pixi-home"
export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
pixi install --locked 2>&1 | tail -25

"$REPO_ROOT/.pixi/envs/default/bin/python" -c \
  "import numpy, scipy; print('numpy', numpy.__version__, '| scipy', scipy.__version__)"
echo INSTALL_DONE