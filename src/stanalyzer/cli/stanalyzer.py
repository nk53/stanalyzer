import argparse
import enum
import io
import os
import re
import sys
from collections.abc import Iterable
from importlib import import_module
from pathlib import Path
from typing import Any, TypeAlias, cast
from .validators import exec_name, p_int, p_float
from ..utils import read_json


def ExistingFile(path: str) -> str:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"No such file: '{path}'")
    if not p.is_file():
        raise ValueError(f"'{path}' is not a regular file")
    return path


class LazyFile:
    _lazyfile_obj: 'FileLike | None'

    def __init__(self, *args, **kwargs):
        """Doesn't open a file until access is attempted. Args eventually
        passed to builtin open()
        """
        self._lazyfile_args = args
        self._lazyfile_kwargs = kwargs
        self._lazyfile_obj = None

    def __enter__(self) -> io.TextIOBase:
        return self._lazyfile_resolve().__enter__()

    def __exit__(self, *args, **kwargs):
        return self._lazyfile_resolve().__exit__(*args, **kwargs)

    def __getattr__(self, name: str) -> Any:
        self._lazyfile_resolve()
        return getattr(self._lazyfile_obj, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith('_lazyfile_'):
            object.__setattr__(self, name, value)

        self._lazyfile_resolve()
        setattr(self._lazyfile_obj, name, value)

    def __iter__(self):
        return iter(self._lazyfile_resolve())

    def _lazyfile_resolve(self) -> io.TextIOBase:
        """Opens the file if not already opened"""
        if self._lazyfile_obj is None:
            # if file cannot be accessed, this may raise an OSError
            self._lazyfile_obj = open(*self._lazyfile_args, **self._lazyfile_kwargs)

        return cast(io.TextIOBase, self._lazyfile_obj)

    def resolve(self) -> io.TextIOBase:
        """Returns the file handle, opening it, if necessary"""
        return self._lazyfile_resolve()


class LazyJSONReader:
    def __init__(self, json_filename):
        self._lazyjsonreader_file = json_filename
        self._lazyjsonreader_data = None

    def __getattr__(self, name):
        self._lazyjsonreader_resolve()
        return getattr(self._lazyjsonreader_data, name)

    def __setattr__(self, name, value):
        if name.startswith('_lazyjsonreader_'):
            return object.__setattr__(self, name, value)

        self._lazyjsonreader_resolve()
        return setattr(self._lazyjsonreader_data, name, value)

    def __getitem__(self, key):
        return self._lazyjsonreader_resolve()[key]

    def __setitem__(self, key, value):
        self._lazyjsonreader_resolve()[key] = value

    def _lazyjsonreader_resolve(self) -> dict:
        if self._lazyjsonreader_data is None:
            self._lazyjsonreader_data = read_json(self._lazyjsonreader_file)
        return cast(dict, self._lazyjsonreader_data)


def add_exec_args(subparser: Any, *args: str) -> None:
    def add_exec_arg(arg):
        default = exec_name(arg, error=False)

        def path_or_raise(value: str) -> str | None:
            if value:
                return exec_name(value)

            if default:
                return default

            raise argparse.ArgumentTypeError(f"No {arg} in PATH")

        if not arg:
            raise ValueError("arg cannot be empty")

        arg_safe = f"--{arg.replace('_ ', '-')}-path"
        if default is None:
            print(f"Warning: no {arg} in PATH", file=sys.stderr)
            subparser.add_argument(arg_safe, metavar='PATH', type=path_or_raise, required=True,
                                   help=f"Path to {arg} executable. "
                                        f"REQUIRED because {arg} is not in PATH.")
        else:
            subparser.add_argument(arg_safe, metavar='PATH', type=path_or_raise,
                                   help=f"Path to {arg} executable. "
                                        f"REQUIRED because {arg} is not in PATH.")

    for arg in args:
        add_exec_arg(arg)


def add_project_args(subparser, *args):
    global from_settings

    defstr = '(default: use project settings)'
    if 'psf' in args:
        subparser.add_argument(
            '-p', '--psf', type=ExistingFile, default=from_settings, metavar='FILE',
            help=f"File containing system topology {defstr}")
    if 'traj' in args:
        subparser.add_argument(
            '-t', '--traj', type=ExistingFile, nargs='+', default=from_settings, metavar='FILE',
            help=f"One or more coordinate containing files {defstr}")
    if 'out' in args:
        subparser.add_argument(
            '-o', '--out', type=argparse.FileType('w'), default=sys.stdout,
            help="File to write results (default: stdout)")
    if 'time_step' in args:
        subparser.add_argument(
            '-ts', '--time_step', type=p_float, default=from_settings,
            help="Amount of time between frames in trajectory files")
    if 'interval' in args:
        subparser.add_argument(
            '-i', '--interval', type=p_int, default=1,
            help="step size when reading frames (default: 1 = read every frame)")


def get_settings(analysis_name: str | None = None) -> dict:
    global analysis, defaults

    parser = argparse.ArgumentParser(
        description='Run or list available analysis modules.', prog='stanalyzer',
        usage="%(prog)s [-hlc] [-s FILE] [analysis_name [-h] [ARGS ...]]")
    if analysis_name is None:
        parser.add_argument('analysis_name', nargs='?', help='Name of analysis to run. '
                            'See `stanalyzer analysis_name -h` for analysis specific options')
    parser.add_argument('-l', '--list', action='store_true',
                        help='List available analysis modules')
    parser.add_argument('-c', '--compact', action='store_true',
                        help='Show available modules on one line. Implies --list.')
    parser.add_argument('-s', '--settings', metavar='file', default='project.json',
                        help='specify name of file containing settings (default: project.json)')

    # some manual parsing is required because parse_known_args() will still
    # fail if an arg is valid both before and after analysis_name
    full_args = sys.argv[1:]

    # get all args before analysis_name; relies on analysis_name being the
    # first positional argument
    if analysis_name is None:
        main_args = []
        prev_arg = None
        idx = 0
        for idx, arg in enumerate(full_args):
            if arg.startswith('-'):
                main_args.append(arg)
            elif prev_arg in ['-s', '--settings']:
                main_args.append(arg)
            else:
                analysis_name = arg
                main_args.append(arg)
                break
            prev_arg = arg

        args = parser.parse_args(main_args)
        remaining_args = full_args[idx+1:] if len(full_args) > idx else []

    if args.list or args.compact:
        print("List of available analysis modules:")
        analyses = sorted(list_analyses())

        if args.compact:
            print(f"    {', '.join(analyses)}")
            sys.exit(0)

        for analysis in analyses:
            print(f"    {analysis}")

        sys.exit(0)

    # TODO: allow only registered imports
    if analysis_name is None:
        analysis_name = args.analysis_name
    if analysis_name is None:
        print("error: missing analysis name")
        parser.print_help()
        sys.exit(15)

    if analysis_name == 'config':
        # run config tool and exit
        from . import config
        config.main()
        sys.exit(0)

    if Path(args.settings).exists():
        defaults = LazyJSONReader(args.settings)
    else:
        print(f"Warning: skipping non-existent settings file '{args.settings}'", file=sys.stderr)

    try:
        _import_name = f"stanalyzer.analysis.{analysis_name}"
        analysis = import_module(_import_name)
    except ModuleNotFoundError as exc:
        if exc.name == _import_name:
            print(f"Error: Not a known analysis type: {analysis_name}", file=sys.stderr)
            parser.print_usage(file=sys.stderr)
            sys.exit(1)

        raise  # analysis_name found, but causes ModuleNotFoundError

    analysis_parser: argparse.ArgumentParser = analysis.get_parser()
    analysis_settings = analysis_parser.parse_args(remaining_args)

    settings = analysis_settings.__dict__
    for name, value in settings.items():
        if value is from_settings:
            try:
                value = defaults[name]
            except KeyError:
                print(f"Missing value for '{name}', which must be specified either on the "
                      f"command line or in settings file '{args.settings}'.", file=sys.stderr)
                analysis_parser.print_usage()
                sys.exit(2)

        if name == 'traj':
            value = get_traj(value)

        settings[name] = value

    return settings


def get_traj(pattern: str | Iterable[str | Path]) -> list[Path]:
    """Given a glob or list of files, returns matching files sorted numerically"""
    re_num = re.compile(r'(\d+)')

    if isinstance(pattern, str):
        pattern = Path().glob(pattern)
    elif isinstance(pattern, list):
        pattern = [Path(p) for p in pattern]

    pattern = cast(Iterable[Path], pattern)

    # get numeric components and compare as list of ints
    return sorted(pattern, key=lambda p: list(int(n) for n in re_num.findall(p.name)))


def resolve_file(file_ref: LazyFile | io.TextIOBase | str, mode: str = 'r') -> io.TextIOBase:
    if isinstance(file_ref, str):
        return cast(io.TextIOBase, open(file_ref, mode))
    if isinstance(file_ref, LazyFile):
        return file_ref.resolve()
    return file_ref


def resolve_ts(time_step: int | float | str) -> float:
    if isinstance(time_step, str):
        return float(time_step.split()[0])
    return float(time_step)


def run_server():
    """Entry point for wrapper scripts"""
    def ip_addr(v: str) -> str:
        parts = v.split('.')
        if len(parts) != 4:
            raise ValueError(f"Invalid IP address: '{v}'")
        for part in parts:
            if (intval := int(part)) < 0 or intval > 255:
                raise ValueError(f"Address component out of range: '{part}'")
        return v

    import stanalyzer
    os.chdir(stanalyzer.__path__[0])

    parser = argparse.ArgumentParser(
        description="Start the ST-Analyzer web server")
    parser.add_argument('-p', '--port', type=int, default=8000,
                        help="(default: 8000)")
    parser.add_argument('--host', type=ip_addr, default='127.0.0.1',
                        help="(default: 127.0.0.1)")
    args = parser.parse_args()

    import uvicorn
    uvicorn.run('stanalyzer.main:app', host=args.host, port=args.port,
                log_level='info')


def list_analyses() -> list[str]:
    """Show module name in ST-Analyzer's analysis directory.

    Modules are listed if they define `main` and `get_parser` functions.
    """
    import stanalyzer
    sta_root = stanalyzer.__path__[0]
    analysis_dir = Path(sta_root) / "analysis"

    analyses: list[str] = []
    for path in analysis_dir.iterdir():
        # skip non-python files
        if not path.name.endswith('.py'):
            continue
        # import may have side effects; just check for 'def main(...'
        has_main = has_parser = False
        with path.open() as file_obj:
            for line in file_obj:
                if line.startswith('def main'):
                    has_main = True
                elif line.startswith('def get_parser'):
                    has_parser = True

        if has_main and has_parser:
            analyses.append(path.stem)

    return analyses


def main(analysis_name: str | None = None, settings: "DictLike | None" = None) -> Any:
    # MDAnalysis and its imports like to throw dev warnings that are hard to disable
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")

    if settings is None:
        settings = get_settings()

    assert analysis is not None, "Failed to import analysis module"

    # yes, this second one really is necessary
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")

    return analysis.main(settings)


# type aliases
DictLike: TypeAlias = LazyJSONReader | dict
FileLike: TypeAlias = LazyFile | io.TextIOBase
FileRef: TypeAlias = FileLike | str
FileRefList: TypeAlias = list[FileRef]

# module containing analysis to run
analysis: Any | None = None
# lazy-loaded settings file containing defaults to use for absent arguments
defaults: DictLike = {}
# signals that a value should be obtained from the default dict
default_action = enum.Enum('default_action', 'from_settings')
from_settings = default_action.from_settings

if __name__ == '__main__':
    main()
