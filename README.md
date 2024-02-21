# Installation

Currently, HPC-GUI assumes you have the following installed and loaded:

1. Conda/Anaconda
2. A conda environment named `balsam`

If not, you should install Balsam by following the instructions shown [here](https://balsam.readthedocs.io/en/latest/user-guide/installation/).

Alternative installation approaches:

## Conda

```bash
# optionally create an environment
conda create --name hpcgui python=3.9
conda activate hpcgui

# install required packages
conda install fastapi==0.78.0 jinja2==3.1.2 uvicorn==0.18.3 pyyaml==6.0 sphinx==6.1.3 sphinx-js
```

Then, edit `.hpcgui_config` and change `conda activate balsam` to `conda activate hpcgui`.

## pip

I have not tested this approach at all so *caveat emptor*.

```bash
pip install fastapi==0.78.0 jinja2==3.1.2 uvicorn==0.18.3 pyyaml==6.0 sphinx==6.1.3
```

Then, edit `run_server.sh` and remove or comment-out this line:
```bash
source .hpcgui_config
```

# Usage

HPC-GUI is in alpha development stage; it is intended to be run and accessed on your personal computer, but can be made accessible over the internet (see below).

To access the GUI via local-only connection:

1. Navigate a terminal to the HPC-GUI directory.
2. Run `./run_server.sh`, which starts a local server session accessible *only* by the bound IP `127.0.0.1`.
3. Open [127.0.0.1:8000](http://127.0.0.1:8000) in your web browser.

To access the GUI via IPv4:

1. Edit `run_server.sh` and add the option `--host 0.0.0.0` to `uvicorn`.
2. Run `./run_server.sh`.
3. Enter your assigned IP or hostname in your browser's address bar, including the port, e.g. `myhostname:8000`.

You will not be able to access the site via hostname unless your network supports hostname broadcasting or your hostname is registered by DNS.

For IPv6 support, use `--host '::'`.

To access without typing `:8000`, change the value of `--port` to `80`.

# Caching Issues

If you are testing recent development changes to static files such as `forms.js` or `style.css`, you browser may be serving you the old version of a file. Some ways you can get around this:
1. Clear your browser's cache of this site
2. Disable caching for this site
3. Open the site in a new private browsing window
4. If using Firefox with developer tools enabled, open the network tab and toggle the "Disable caching" option.
