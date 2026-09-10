#!/usr/bin/env bash
# One-time setup: install pixi and build the locked linux-64 environment.
#
# Intended to run INSIDE a docker container with a copy of the repo mounted
# somewhere (see README.md, "Dockerized linux-64 testing"). The script derives
# the repo root from its own location, so the mount point does not matter.
set -euxo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # <root>/src/stanalyzer/tests

export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq curl ca-certificates build-essential >/dev/null

curl -fsSL -o /tmp/pixi.tar.gz \
  https://github.com/prefix-dev/pixi/releases/download/v0.73.0/pixi-x86_64-unknown-linux-musl.tar.gz
tar -xzf /tmp/pixi.tar.gz -C /usr/local/bin
pixi --version

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