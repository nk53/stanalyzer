# Installation

Currently, ST-Analyzer should be installed through Conda/Anaconda. Conda's installation instructions begin [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

ST-Analyzer has only been tested on macOS 14.5. It is unlikely to work on Windows in its current state, but may work on Linux.

## Conda installation

You should install all ST-Analyzer dependencies in a separate environment so that they do not conflict with any other versions of packages you have. All relevant dependencies are included in `envs/sta-dev.yml`. To install them in a new `sta-dev` environment, run this command:

```bash
conda env create -f envs/sta-dev.yml
```

# Usage

ST-Analyzer is in alpha development stage; it is intended to be run and accessed on your personal computer, but can be made accessible over the internet (see below).

## The long way

This approach should work without any hiccups. If it does not, please refer back to the installation instructions and ensure that you installed ST-Analyzer correctly.

First, navigate your terminal to the ST-Analyzer directory. Then, run these commands:

```bash
conda activate sta-dev
uvicorn main:app --port 8000 --reload
```

This will run a web server in your shell's foreground and show you several debugging messages. The server will listen for connections on port 8000 and also automatically restart if any Python or Jinja2 files change while it is running. Note that it does not detect changes in other files (YAML, JS, CSS, etc.); to reflect changes in those files, you will need to manually restart the server.

The server can be terminated by sending a keyboard interrupt signal (Ctrl-C on Mac/Linux) to its controlling terminal. To start it again, simply rerun the `uvicorn` command.

Access the server by opening [127.0.0.1:8000](http://127.0.0.1:8000) in your web browser.

## The short way (read long way first)

If you're like me, you don't want to type more than five keys to do a common task like starting the server (especially if you need to restart it frequently) if you don't have to. The script `./run_server.sh` checks if the sta-dev environment is already loaded, loads it if necessary, starts the server, and (if using macOS) copies the web address+port to your clipboard.

It has not been tested extensively, so YMMV. Check if it works like so (from ST-Analyzer root directory):

```bash
./run_server.sh
```

Since no other file in this directory should start with an `r`, this can be accomplished with five key presses: `.`, `/`, `r`, `<Tab>`, `<Enter>`. If everything worked correctly, you should be able to paste the address into your browser's address bar.

## Allowing remote connections

To access the GUI via IPv4:

1. Edit `run_server.sh` and add the option `--host 0.0.0.0` to `uvicorn`.
2. Run `./run_server.sh`.
3. Enter your assigned IP in your browser's address bar, including the port, e.g. `my.local.ip.addr:8000`.

You will not be able to access the site via hostname unless your network supports hostname broadcasting or your hostname is registered by DNS.

For IPv6 support, use `--host '::'`.

To access without typing `:8000`, change the value of `--port` to `80`.

# Caching Issues

If you are testing recent development changes to static files such as `forms.js` or `style.css`, you browser may be serving you the old version of a file. Some ways you can get around this:
1. Clear your browser's cache of this site
2. Disable caching for this site
3. Open the site in a new private browsing window
4. If using Firefox with developer tools enabled, open the network tab and toggle the "Disable caching" option.
