#!/bin/bash
set -e

# unset previous conda vars
unset CONDA_SHLVL
unset _CE_CONDA
unset _CE_M
unset CONDA_EXE

# does this shell have conda in its PATH?
echo "Trying to find conda ..."
if ! installer=$(which conda); then
    # does the default python have a conda module?
    if python -m conda -h; then
        installer="python -m conda"
    fi
fi >/dev/null 2>/dev/null

if [ -z "$installer" ]; then
    echo "Can't find conda executable."
    echo "Please enter path to conda executable, or leave the line blank and hit Enter to cancel."
    prompt="Path to conda executable (default: quit): "
    while [ 1 ]; do
        read -p "$prompt" installer
        [ -z "$installer" ] && exit 1
        if [ ! -e "$installer" ]; then
            echo "Sorry, '$installer' does not exist. Try again or hit Enter to cancel."
        elif [ ! -f "$installer" ]; then
            echo "Sorry, '$installer' is not a regular file. Try again or hit Enter to cancel."
        elif [ ! -x "$installer" ]; then
            echo "Sorry, '$installer' is not executable. Try again or hit Enter to cancel."
        else
            # got an actual program, at least
            break
        fi
    done
fi

# initialize conda
echo "found '$installer'. Trying to initialize ..."
eval "$($installer 'shell.bash' 'hook')"

# ensure conda and conda-build are up-to-date
echo "Checking that conda and conda-build are up to date ..."
conda install -y conda-build
conda update -y conda conda-build

# download repo and run build
orig_dir=sta-dir
clone_dir=$orig_dir
if [ -e "$clone_dir" ]; then
    echo "Can't clone into '$clone_dir' because it already exists."
    echo "Please remove '$clone_dir' or enter a new path below."
    echo "Enter a blank line to cancel."
    while [ 1 ]; do
        read -p "Path to save stanalyzer source: " clone_dir
        [ -z "$clone_dir" ] && exit 1
        if [ "$clone_dir" == "$orig_dir" ]; then
            echo "Please remove '$orig_dir' and try again."
            exit 1
        elif [ -e "$clone_dir" ]; then
            echo "$clone dir already exists. Please try again or hit Enter to cancel."
        else
            break
        fi
    done
fi

git clone "https://github.com/nk53/stanalyzer" "$clone_dir"
cd "$clone_dir"/conda.build
./local_make.sh
cd ../..

# remove repo
chmod -R 777 "$clone_dir"
rm -r "$clone_dir"
