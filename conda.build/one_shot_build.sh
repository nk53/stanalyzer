#!/bin/bash
set -e
conda init
conda activate base
conda install -y conda-build
conda update -y conda conda-build

git clone "https://github.com/nk53/stanalyzer" sta-build
cd sta-build/conda.build
./local_make.sh
cd ../..
rm -r sta-build
