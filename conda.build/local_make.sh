#!/bin/bash

output_dir=stanalyzer-builds
output_dir_parent="$HOME/Downloads"
output_abspath="$output_dir_parent/$output_dir"
output_archive="$output_dir.tgz"

if [ ! -e "$output_abspath" ]; then
    mkdir -p "$output_abspath"
fi

CHANNEL="-c conda-forge"
OPTS="--keep-going --no-anaconda-upload"

set -e

conda build --variants '{suffix: ["", client]}' \
    --output-folder "$output_abspath" $OPTS $CHANNEL .
conda build --no-remove-work-dir --variants "{suffix: [dev]}" \
    --output-folder "$output_abspath" $OPTS $CHANNEL .

echo "cleaning up intermediate files ..."
conda build purge

echo "build complete, compressing ..."
cd "$output_dir_parent"
tar cvzf "$output_archive" "$output_dir"
rm -r "$output_dir"

echo "final output: $output_dir_parent/$output_archive"
