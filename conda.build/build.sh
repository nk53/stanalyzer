if [ "$STA_SUFFIX" == "dev" ]; then
    # just copy files necessary for post-installer
    mkdir $SP_DIR/stanalyzer
    cp $SRC_DIR/conda.build/install.py $SP_DIR/stanalyzer
    cp $SRC_DIR/conda.build/__init__.py $SP_DIR/stanalyzer

else
    if [ -z "$target_platform" ]; then
        echo missing target platform
        exit 1
    fi
    python -m pip install --no-deps --ignore-installed -vv .
fi

pth="$HOME/local/site-packages/stanalyzer"
if [ -h "$pth" ]; then
    echo "WARNING: previous stanalyzer installation conflicts with this one"
    echo "Moving old installation:"
    mv -v "$pth" "$pth.bak"
fi

