#!/bin/bash
set -e

# unset previous conda vars
unset CONDA_SHLVL
unset _CE_CONDA
unset _CE_M
unset CONDA_EXE

# initialize conda
eval "$(python -m conda shell.bash hook)"

# ensure conda and conda-build are up-to-date
conda activate base
conda install -y conda-build
conda update -y conda conda-build

# download repo and run build
git clone "https://github.com/nk53/stanalyzer" sta-build
cd sta-build/conda.build
./local_make.sh
cd ../..

# remove repo
rm -r sta-build
