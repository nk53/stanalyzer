if [ "$STA_SUFFIX" == "dev" ]; then
    python -m pip install --no-deps --ignore-installed --editable \
        "git+https://github.com/nk53/stanalyzer.git@sta#egg=stanalyzer"
else
    python -m pip install --no-deps --ignore-installed .
fi

pth="$HOME/local/site-packages/stanalyzer"
if [ -h "$pth" ]; then
    echo "WARNING: previous stanalyzer installation conflicts with this one"
    echo "Moving old installation:"
    mv -v "$pth" "$pth.bak"
fi

