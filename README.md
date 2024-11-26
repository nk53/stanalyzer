# Installation

Currently, ST-Analyzer should be installed through Conda/Anaconda. Conda's installation instructions begin [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

ST-Analyzer has been tested on macOS (x64 and arm64), Ubuntu (x64), and WSL (x64, running Ubuntu).

See the section below for regular installation.

For development installation, skip to the section __Installation via Conda (ST-Analyzer developers)__.

## Installation via Conda (regular users)

ST-Analyzer works best when it is installed in a fresh conda environment using the `conda-forge` channel. The example below installs into a new environment named `sta`:

```bash
conda create -n sta -c conda-forge nk53::stanalyzer
```

To use ST-Analyzer, first activate its environment:

```bash
conda activate sta  # required once per shell session
```

Then, you should have two new commands available: `stanalyzer` and `sta-server`. Their usage is explained in the __Usage__ section.

### Troubleshooting: Installation
If you get a message like "environment already exists", you should either use a different name or remove the old one. This command removes the sta environment:

```bash
# WARNING: Removes ST-Analyzer's environment (or other w/ name sta)
conda env remove -n sta
```

## Installation via Conda (ST-Analyzer developers)

Development installation entails obtaining the ST-Analyzer source code and linking it to your python environment. Fortunately, this is easy to do within conda's framework:

```bash
# install ST-Analyzer's dependencies and post-install helper
conda create -n sta-dev -c conda-forge nk53::stanalyzer-dev

# enter environment
conda activate sta-dev

# download source code and link it in the environment's site-packages
stanalyzer-finish-install
```

Without any options, `stanalyzer-finish-install` asks where to save the source code before downloading and linking anything.

The helper script is just an interactive version of these two commands:

```bash
git clone "git@github:nk53/stanalyzer.git" dir_to_save_to
pip install --no-deps --ignore-installed --editable dir_to_save_to
```

This is known as an "editable" installation; updates to the source code reflect immediately without the need to repackage/reinstall `stanalyzer-dev`. The only reason to mess with conda is if something in the upstream environment changes, in which case `conda update stanalyzer-dev` should be enough.

### Caching Issues

If you are testing recent development changes to static files such as `forms.js` or `style.css`, you browser may be serving you the old version of a file. Some ways you can get around this:
1. Clear your browser's cache of this site
2. Disable caching for this site
3. Open the site in a new private browsing window
4. If using Firefox with developer tools enabled, open the network tab and toggle the "Disable caching" option.

# Usage

ST-Analyzer is in alpha development stage; it is intended to be run and accessed on your personal computer.

## Web GUI: sta-server

If ST-Analyzer was installed successfully, you should be able to start the Web GUI with 1-2 commands, run from anywhere:

```bash
conda activate sta  # if not already active
sta-server
```

This runs the server process in the foreground of your current shell.

If you did not receive any error messages, then the server should be listening for localhost connections, and will show the corresponding URL. At the time of this writing, it should be `http://127.0.0.1:8000`. You can copy and paste that into your browser's address bar to reach the ST-Analyzer front page.

## Command-line tool: stanalyzer

Advanced users often want to automate analysis themselves without dealing with the web UI. ST-Analyzer was designed with this goal in mind. All analysis processes can be run from the command line directly. In fact, this is what the web GUI is doing behind the scenes.

The easiest way to get started with the command-line interface (CLI) is to first use the web interface to create a project and run a small sample analysis, e.g., the "System size" analysis with default settings. When an analysis job is successfully submitted through the web GUI, it writes the exact command needed to run the analysis from the command-line. E.g.:

```bash
args: stanalyzer system_size --time_step 1 --out system_size.out
```

In all cases, the web GUI uses the verbose version of a command's options. To see all available options for a command and their aliases, use `stanalyzer [analysis_program] -h`. E.g.:

```bash
stanalyzer system_size -h
```
