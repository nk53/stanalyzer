#!/bin/bash

output_folder=$HOME/Downloads/stanalyzer-builds
output_archive=$output_folder.tgz

if [ ! -e "$output_folder" ]; then
    mkdir "$output_folder"
fi

CHANNEL="-c conda-forge"
if [ "$1" == "-f" ]; then
    OPTS="--keep-going"
else
    OPTS="--skip-existing"
fi
OPTS="$OPTS --output-folder \"$output_folder\" --no-anaconda-upload"

set -e

conda build --variants '{suffix: ["", client]}' $OPTS $CHANNEL .
conda build --no-remove-work-dir --variants "{suffix: [dev]}" $OPTS $CHANNEL .

echo "cleaning up intermediate files ..."
conda build purge

echo "build complete, compressing ..."
tar cvzf "$output_archive" "$output_folder"
rm -r "$output_folder"

echo "final output: $output_archive"
