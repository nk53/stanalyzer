python -m pip install --no-deps --ignore-installed .

pth="$HOME/local/site-packages/stanalyzer"
if [ -h "$pth" ]; then
    echo "WARNING: previous stanalyzer installation conflicts with this one"
    echo "Moving old installation:"
    mv -v "$pth" "$pth.bak"
fi

