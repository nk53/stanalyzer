#!/bin/bash

if grep -q 'keep-egg' <<< $@; then
    keep=1
fi

export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
pip install -e .
echo
