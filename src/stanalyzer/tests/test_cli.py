"""Compare results vs. previous runs"""
import io
import re
import typing as t
import unittest
from collections.abc import Callable
from pathlib import Path

import invoke

from stanalyzer.validation import Project
from stanalyzer.utils import write_settings

T = t.TypeVar('T')
P = t.ParamSpec('P')
TestFunction: t.TypeAlias = Callable[[Callable[P, T]], Callable[P, T]]
StrToPath: t.TypeAlias = Callable[[str | Path], Path]
IOPair: t.TypeAlias = tuple[io.TextIOWrapper, io.TextIOWrapper]
IOTriple: t.TypeAlias = tuple[io.TextIOWrapper, io.TextIOWrapper, str]
IO2_3: t.TypeAlias = IOPair | IOTriple

capsword_start = re.compile(r'([^a-z])')


def camel_to_snake(name: str) -> str:
    """Takes NameLikeThis and returns name_like_this, if possible"""
    if not (match := capsword_start.split(name)):
        return name

    # skip possible leading ''
    while match and match[0] == '':
        match.pop(0)

    words: list[str] = []
    while match:
        start = match.pop(0)
        rest = match.pop(0)
        words.append(start + rest)

    return '_'.join(words).lower()


def skipUnlessAttrNotNone(obj: object, attr: str) -> TestFunction:
    if getattr(obj, attr, None) is None:
        return unittest.skip(f"{obj!r} doesn't have {attr!r}")
    return lambda func: func


class ManagedConfig:
    """Dynamically creates project.json configs.

    If used as a context manager, the file is removed on context exit.
    """

    def __init__(self, *,
                 title: str = "Test Case",
                 input_relpath: str | Path,
                 output_relpath: str | Path,
                 traj: str,
                 psf: str,
                 time_step: str = "1 ns",
                 scheduler: str = "interactive"):
        inp = Path('.').resolve() / input_relpath
        out = Path('.').resolve() / output_relpath

        project = Project(
            title=title,
            input_path=inp,     # absolute path
            output_path=out,    # absolute path
            traj=traj,          # relative to input_path
            psf=psf,            # relative to input_path
            time_step=time_step,
            scheduler=scheduler,
        )

        project_dict = project.model_dump(mode="json")
        project_dict.pop('id', None)

        # Project does str -> dict for validation; revert to str for writing
        project_dict['time_step'] = time_step

        self.inp = inp
        self.out = out
        self.config = project
        self.config_path = out / 'project.json'
        self.project_dict = project_dict

    def write(self) -> None:
        """Write project.json to output dir.

        Overwrites existing file, if present.
        """
        write_settings(path=self.out / 'project.json', data=self.project_dict)

    def __enter__(self) -> Project:
        self.write()
        return self.config

    def __exit__(self, x: t.Any, y: t.Any, z: t.Any) -> t.Literal[False]:
        self.config_path.unlink()
        return False


class AnalysisCase(unittest.TestCase):
    ctx: invoke.Context
    config_path: Path
    config: Project
    manager: ManagedConfig

    def __init__(self, methodName='runTest'):
        super().__init__(methodName)

        manager: ManagedConfig | None = getattr(self, 'manager', None)
        if manager is None:
            return

        self.config = manager.config

        outdir = Path(self.config.output_path)

        if outdir.is_file():
            raise FileExistsError(f"{outdir} exists and is not a directory")

        if not outdir.exists():
            outdir.mkdir(parents=True)

    def infile(self, relpath: str | Path) -> Path:
        return self.config.input_path / relpath

    def outfile(self, relpath: str | Path) -> Path:
        return self.config.output_path / relpath

    def file_empty(self, path: str | Path | io.TextIOWrapper,
                   path_type: StrToPath | None = None) -> bool:
        return self.file_size(path, path_type) == 0

    def file_exists(self, path: str | Path, path_type: StrToPath | None = None) -> bool:
        if path_type is None:
            return Path(path).exists()
        return path_type(path).exists()

    def file_size(cls, path: str | Path | io.TextIOWrapper,
                  path_type: StrToPath | None = None) -> int:
        if isinstance(path, (str, Path)):
            abspath = Path(path)
            if not abspath.is_absolute():
                if path_type:
                    abspath = path_type(abspath)

            with open(abspath) as outfile:
                outfile.seek(0, io.SEEK_END)
                return outfile.tell()
        else:
            outfile = path
            outfile.seek(0, io.SEEK_END)
            return outfile.tell()

    @t.overload
    def run_analysis(self, args: str, stem: str = '', analysis: str = '',
                     accepts_o: t.Literal[True] = True) -> IOTriple: ...

    @t.overload
    def run_analysis(self, args: str, stem: str = '', analysis: str = '',
                     accepts_o: t.Literal[False] = False) -> IOPair: ...

    @t.overload
    def run_analysis(self, args: str, stem: str = '', analysis: str = '',
                     accepts_o: bool = ...) -> IO2_3: ...

    def run_analysis(self, args: str, stem: str = '', analysis: str = '',
                     accepts_o: bool = True) -> IO2_3:
        """Shortcut for when out and err share a stem with analysis name

        If `stem` is empty, self.analysis_name is used.
        If `analysis` is empty, then it is assumed to be the same as `stem`.
        If `accepts_o`, then --out {name}.dat is added to args.
        """

        if not stem:
            assert hasattr(self, 'analysis_name'), \
                "Can't omit both stem and self.analysis_name"
            stem = self.analysis_name

        if not analysis:
            analysis = getattr(self, 'analysis_name', stem)

        if accepts_o:
            datfile = f'{stem}.dat'
            dat_abspath = str(self.outfile(datfile))
            args = ' '.join([f'--out "{datfile}"', args])

        out_err = self.run_cmd(f'stanalyzer {analysis} {args}',
                               out_filename=f'{stem}.out',
                               err_filename=f'{stem}.err')

        if accepts_o:
            return out_err + (dat_abspath,)
        return out_err

    def run_cmd(self, args: str, out_filename: str, err_filename: str) -> IOPair:
        out_path = self.config.output_path / out_filename
        err_path = self.config.output_path / err_filename

        out_stream = out_path.open('w')
        err_stream = err_path.open('w')
        print('args:', args, file=out_stream)
        out_stream.flush()

        self.ctx.run(args, out_stream=out_stream, err_stream=err_stream)

        return out_stream, err_stream

    def run(self, result: unittest.TestResult | None = None) -> unittest.TestResult | None:
        if not getattr(self, '__unittest_skip__', False):
            self.ctx = invoke.Context()

            assert self.manager is not None, "Missing project.json"

            self.manager.write()
            with self.ctx.cd(self.config.output_path):
                result = super().run(result)

        return result


class SoohyungCase(AnalysisCase):
    """Shortcut for preparing project.json using soohyung_membrane as the template"""
    default_output: t.ClassVar[str | Path] = 'test_case'
    standard_args: t.ClassVar[str | None] = None
    test_standard: Callable

    # subclass should override if its name doesn't follow camel_to_snake scheme
    analysis_name: t.ClassVar[str] = ''

    def standard_test(self) -> None:
        outfile = self.outfile
        args = self.standard_args

        assert args is not None

        out, err, dat = self.run_analysis(args)

        self.assertTrue(self.file_exists(dat, outfile))
        self.assertFalse(self.file_empty(out, outfile))

    def __init_subclass__(cls, **kwargs):
        if not cls.analysis_name:
            cls.analysis_name = camel_to_snake(cls.__name__)
        cls.default_output = Path('results') / cls.analysis_name
        cls.test_standard = skipUnlessAttrNotNone(cls, 'standard_args')(SoohyungCase.standard_test)

        super().__init_subclass__(**kwargs)

    def __init__(self, methodName='runTest'):
        if not hasattr(self, 'manager'):
            self.manager = ManagedConfig(
                input_relpath="soohyung_membrane",
                output_relpath=self.default_output, traj="step7_*.dcd",
                psf="step5_input.psf")
        super().__init__(methodName)


class DensityZ(SoohyungCase):
    standard_args = '-c --sel "name C*"'


class ContactResidenceTime(SoohyungCase):
    standard_args = '--sel "protein and name CA" --threshold "5.0"'
    analysis_name = 'contact_res_time'


class HBond(SoohyungCase):
    standard_args = '--sel "segid MEMB or protein" --hydrogens-sel "None" ' \
           '--acceptors-sel "None" --d-a-cutoff "3.0" ' \
           ' --d-h-a-angle-cutoff "150.0"'
    analysis_name = 'hbond'


class PiStacking(SoohyungCase):
    standard_args = '--sel "not segid TIP3 and not segid IONS" ' \
           ' --pi-pi-dist-cutoff "6.0" --pi-cation-dist-cutoff "6.0"'


class RMSF(SoohyungCase):
    analysis_name = 'rmsf'
    sel = "segid PROA and name CA"
    standard_args = f'--sel-align "{sel}" --sel-rmsf "{sel}"'


@unittest.skip("missing test files & settings")
class SaltBridge(AnalysisCase):
    @unittest.skip("missing test files & settings")
    def test_standard(self) -> None:
        outfile = self.outfile

        sel = "not segid TIP3 and not segid IONS"
        pos = "resname ARG LYS and name NE NH* NZ*"
        neg = "resname ASP GLU and name OE* OD*"
        args = f'--sel "{sel}" --negative-sel "{sel}" --positive-def "{pos}"'\
               f'--negative-def "{neg}" --dist-cutoff "4.5"'

        out, err, dat = self.run_analysis(args)
        self.assertTrue(self.file_exists(dat, outfile))
        self.assertFalse(self.file_empty(out, outfile))


class SystemSize(SoohyungCase):
    standard_args = ''


class Thickness(SoohyungCase):
    standard_args = '-c --sel "segid MEMB and (name P or name N or name C1[0-9] or name O1[0-9])"'


class WaterBridge(SoohyungCase):
    standard_args = '--sel "protein" --sel2 "None" --water-sel "resname TIP3" '\
                    '--d-a-cutoff "3.0" --d-h-a-angle-cutoff "150.0"'


if __name__ == '__main__':
    unittest.main()
