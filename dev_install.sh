#!/bin/bash

if grep -q 'keep-egg' <<< $@; then
    keep=1
fi

export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0dev0
pip install -e .
echo

if [ ! $keep ]; then
    egg_dir=stanalyzer.egg-info
    [ -d "$egg_dir" ] && rm -r "$egg_dir"
fi

echo "To uninstall development version, run this command:"
echo "    $ pip uninstall stanalyzer"
