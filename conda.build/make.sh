#!/bin/bash

CHANNEL="-c conda-forge"
if [ "$1" == "-f" ]; then
    OPTS="--keep-going"
else
    OPTS="--skip-existing"
fi

set -e

conda build --variants '{suffix: ["", client]}' $OPTS $CHANNEL .
conda build --no-remove-work-dir --variants "{suffix: [dev]}" $OPTS $CHANNEL .
conda build purge
