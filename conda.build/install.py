import os
import sys

pip_e_URL = "https://pip.pypa.io/en/stable/topics/local-project-installs" +\
            "/#editable-installs"
repo = "git@github.com:nk53/stanalyzer.git"


def finish_install():
    saveas = 'stanalyzer'

    prog, args = os.path.basename(sys.argv.pop(0)), sys.argv
    match args:
        case []:
            print("This will clone the git repository and install it as an 'editable' package.")
            print("Read here for more info:")
            print("   ", pip_e_URL)
            print()
            try:
                saveas = input(f"Desired repository save path (default: {saveas}): ") or saveas
            except KeyboardInterrupt:
                print()
                sys.exit(1)
        case '-q':
            pass
        case '-h' | '--help':
            print(f"Usage: {prog} [-h] [-q | -p SAVE_PATH]", file=sys.stderr)
            sys.exit(0)
        case ['-p', saveas]:
            pass
        case _:
            print(f"Usage: {prog} [-h] [-q | -p SAVE_PATH]", file=sys.stderr)
            sys.exit(1)

    cmds = [
        f"git clone {repo} {saveas}",
        f"pip install --no-deps --ignore-installed --editable {saveas}",
    ]

    for cmd in cmds:
        print("Running:", cmd)
        result = os.system(cmd)
        if result:
            print(cmd[:10], '... failed with status', result, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(result)

    print("Done")
