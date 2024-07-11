import argparse
import enum
import io
import re
import sys
from collections.abc import Iterable
from importlib import import_module
from pathlib import Path
from typing import Any, Optional, TypeAlias, cast

from ..utils import read_json


class ExistingFile:
    def __new__(cls, path: str):
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"No such file: '{path}'")
        if not p.is_file():
            raise ValueError(f"'{path}' is not a regular file")
        return path


class LazyFile:
    def __init__(self, *args, **kwargs):
        """Doesn't open a file until access is attempted. Args eventually
        passed to builtin open()
        """
        self._lazyfile_args = args
        self._lazyfile_kwargs = kwargs
        self._lazyfile_obj = None

    def __enter__(self):
        self._lazyfile_resolve()
        return self._lazyfile_obj.__enter__()

    def __exit__(self, *args, **kwargs):
        self._lazyfile_resolve()
        return self._lazyfile_obj.__exit__(*args, **kwargs)

    def __getattr__(self, name):
        self._lazyfile_resolve()
        return getattr(self._lazyfile_obj, name)

    def __setattr__(self, name, value):
        if name.startswith('_lazyfile_'):
            return object.__setattr__(self, name, value)

        self._lazyfile_resolve()
        return setattr(self._lazyfile_obj, name, value)

    def __iter__(self):
        self._lazyfile_resolve()
        return iter(self._lazyfile_obj)

    def _lazyfile_resolve(self):
        """Opens the file if not already opened"""
        if self._lazyfile_obj is None:
            # if file cannot be accessed, this may raise an OSError
            self._lazyfile_obj = open(*self._lazyfile_args, **self._lazyfile_kwargs)

    def resolve(self) -> io.TextIOWrapper:
        """Returns the file handle, opening it, if necessary"""
        self._lazyfile_resolve()
        return self._lazyfile_obj


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
        self._lazyjsonreader_resolve()
        return self._lazyjsonreader_data[key]

    def __setitem__(self, key, value):
        self._lazyjsonreader_resolve()
        self._lazyjsonreader_data[key] = value

    def _lazyjsonreader_resolve(self):
        if self._lazyjsonreader_data is None:
            self._lazyjsonreader_data = read_json(self._lazyjsonreader_file)


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


def get_settings(analysis_name: Optional[str] = None) -> dict:
    global analysis, defaults

    parser = argparse.ArgumentParser(description='TODO', prog='stanalyzer',
                                     add_help=False)
    if analysis_name is None:
        parser.add_argument('analysis_name')
    parser.add_argument('-s', '--settings', metavar='file', default='project.json',
                        help='specify name of file containing settings (default: project.json)')

    # show this parser's help ONLY if it is the first arg
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
        parser.print_help(file=sys.stderr)
        sys.exit()

    args, remaining_args = parser.parse_known_args()

    if Path(args.settings).exists():
        defaults = LazyJSONReader(args.settings)
    else:
        print(f"Warning: skipping non-existent settings file '{args.settings}'", file=sys.stderr)

    # TODO: allow only registered imports
    if analysis_name is None:
        analysis_name = args.analysis_name

    try:
        analysis = import_module(f"stanalyzer.analysis.{analysis_name}")
    except ModuleNotFoundError:
        # TODO: show registered analysis modules in a convenient list
        print(f"Error: Not a known analysis type: {analysis_name}", file=sys.stderr)
        parser.print_usage(file=sys.stderr)
        sys.exit(1)

    analysis_parser = analysis.get_parser()
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


def main(analysis_name: Optional[str] = None, settings: Optional["DictLike"] = None) -> Any:
    if settings is None:
        settings = get_settings()

    assert analysis is not None, "Failed to import analysis module"

    # MDAnalysis likes to throw out warnings that are hard to disable when irrelevant
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
analysis = None
# lazy-loaded settings file containing defaults to use for absent arguments
defaults: DictLike = {}
# signals that a value should be obtained from the default dict
from_settings = enum.Enum('default_action', 'from_settings')

if __name__ == '__main__':
    main()
