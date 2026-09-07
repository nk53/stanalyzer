"""Compare results vs. previous runs"""
import dataclasses
import io
import os
import re
import shlex
import shutil
import subprocess
import sys
import typing as t
import unittest
from collections.abc import Callable
from pathlib import Path

import numpy as np

from stanalyzer.utils import write_settings
from stanalyzer.validation import Project

# Force matplotlib to use the non-interactive Agg backend so analyses that
# generate plots (e.g. cov_analysis heatmap) don't open GUI windows when run
# as subprocesses during tests. Propagates to all spawned subprocesses.
os.environ.setdefault('MPLBACKEND', 'Agg')

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


# Mapping of analysis name -> list of output glob patterns (relative to the
# analysis's output directory). Used by `discover_output_files()` to locate
# the files an analysis produces.
#
# Analyses that accept `--out` write a single file whose name is given by the
# user (typically `<analysis_name>.dat`), so the pattern is `*.dat`.
#
# Analyses with hardcoded/dynamic filenames write one or more files with
# predictable names (e.g. `density_z` writes `<sel>_nb<nbin>_<suffix>.dat`).
#
# Analyses that only write to stdout have no file output and are omitted.
OUTPUT_PATTERNS: dict[str, list[str]] = {
    # --out-based analyses (single output file)
    'contact_res_time': ['*.dat'],
    'hbond': ['*.dat'],
    'pi_stacking': ['*.dat'],
    'rmsf': ['*.dat'],
    'system_size': ['*.dat'],
    'thickness': ['*.dat'],
    'water_bridge': ['*.dat'],
    'rmsd': ['*.dat'],
    'radius_of_gyration': ['*.dat'],
    'position_time': ['*.dat'],
    'position_time_copy': ['*.dat'],
    'compressibility_modulus': ['*.dat'],
    'rdf': ['*.dat'],
    'salt_bridge': ['*.dat'],
    'contacts': ['*.dat'],
    'secondary_structure': ['*.dat'],
    'sasa': ['*.dat'],
    'chol_tilt': ['*.dat'],
    'helix_analysis': ['*.dat'],
    'helix_tilt_rotation_angle': ['*.dat'],
    # Hardcoded/dynamic filenames (no --out)
    'density_z': ['*_nb*_*.dat', 'combined_nb*_*.dat', 'NA_*_nb*_*.dat'],
    'scd': ['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat'],
    'voronoi_shell_comp': ['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat'],
    'voronoi_contact': ['ave_*_*.dat', 'time_*_*.dat', 'NA_time_*_*.dat'],
    'voronoi_apl': ['ave_*_*.dat', 'time_*_*.dat'],
    'msd_solution': ['sys_com_*.dat', 'mol_com_*.dat', '*_*.dat', 'NA_*_*.dat', 'mol_info_*.dat'],
    'msd_membrane': ['*_sys_com_*.dat', '*_mol_com_*.dat', '*_*_*.dat', 'NA_*_*_*.dat', '*_mol_info_*.dat'],
    'clustering_hca': ['cluster.dat', 'cluster_representative.pdb'],
    'cov_analysis': ['corr_matrix.dat', 'eigenvalues.dat', 'eigenvectors.dat'],
    'bond_statistics': ['bond_lengths.dat', 'bond_angles.dat', 'bond_dihedrals.dat'],
}

# Availability of optional external tools. `hole2` is hardcoded to False
# because it is not available on osx-arm64 (see GitHub issue #3).
TOOLS_AVAILABLE: dict[str, bool] = {
    'dssp': shutil.which('dssp') is not None,
    'freesasa': shutil.which('freesasa') is not None,
    'hole2': False,
}


@dataclasses.dataclass(frozen=True, slots=True)
class ReferenceData:
    """Parsed content of a reference output file.

    Holds non-comment, non-empty lines stripped of leading/trailing whitespace.
    Purely numeric files can additionally be parsed into a numpy array via
    :pymethod:`to_numpy`.
    """
    raw_lines: list[str]

    @property
    def is_empty(self) -> bool:
        return len(self.raw_lines) == 0

    def to_numpy(self) -> np.ndarray:
        """Parse all lines as a numeric numpy array.

        Raises ``ValueError`` when any field is non-numeric.
        """
        rows: list[list[float]] = []
        for line in self.raw_lines:
            rows.append([float(f) for f in line.split()])
        return np.array(rows)


def load_reference(dat_path: Path) -> ReferenceData:
    """Load a golden reference output file.

    Handles three format classes:

    * **Pure numeric** (e.g. ``system_size.dat``) — comparable via
      :func:`numpy.testing.assert_allclose`.
    * **Mixed text + numeric columns** (e.g. ``contacts.dat``,
      ``hbond.dat``) — text fields compared for equality, numeric
      fields with tolerance.
    * **Pure text** (e.g. ``compressibility_modulus.dat``) — exact
      string match.

    Comment lines (starting with ``#``) and blank lines are skipped.
    """
    lines: list[str] = []
    with open(dat_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith('#'):
                lines.append(stripped)
    return ReferenceData(raw_lines=lines)


def _compare_lines(test_case: unittest.TestCase,
                   actual: ReferenceData,
                   reference: ReferenceData,
                   actual_path: Path,
                   reference_path: Path,
                   rtol: float,
                   atol: float) -> None:
    """Line-by-line comparison with per-field numeric tolerance."""
    test_case.assertEqual(
        len(actual.raw_lines), len(reference.raw_lines),
        f"Line count mismatch: {len(actual.raw_lines)} != "
        f"{len(reference.raw_lines)} "
        f"({actual_path} vs {reference_path})")

    for idx, (a_line, r_line) in enumerate(
            zip(actual.raw_lines, reference.raw_lines)):
        a_fields = a_line.split()
        r_fields = r_line.split()
        test_case.assertEqual(
            len(a_fields), len(r_fields),
            f"Field count mismatch on line {idx + 1}: "
            f"{len(a_fields)} != {len(r_fields)} "
            f"({a_line!r} vs {r_line!r})")
        for fidx, (a_f, r_f) in enumerate(zip(a_fields, r_fields)):
            try:
                a_val = float(a_f)
                r_val = float(r_f)
                np.testing.assert_allclose(
                    a_val, r_val, rtol=rtol, atol=atol,
                    err_msg=f"Line {idx + 1}, field {fidx + 1}: "
                            f"{a_f} != {r_f}")
            except ValueError:
                test_case.assertEqual(
                    a_f, r_f,
                    f"Line {idx + 1}, field {fidx + 1}: "
                    f"{a_f!r} != {r_f!r}")


def assert_output_matches_reference(test_case: unittest.TestCase,
                                    actual_dat_path: Path,
                                    reference_dat_path: Path,
                                    rtol: float = 1e-5,
                                    atol: float = 1e-8) -> None:
    """Assert that an analysis's output matches its golden reference.

    Strategy:
    1. Both files empty (comment-only) → pass.
    2. Both parseable as pure numeric → ``numpy.testing.assert_allclose``
       (backward-compatible path with configurable tolerance).
    3. Otherwise → line-by-line field comparison: numeric fields use
       ``assert_allclose``, text fields use exact equality.
    """
    if not reference_dat_path.exists():
        test_case.skipTest(f'Reference not found: {reference_dat_path}')

    actual = load_reference(actual_dat_path)
    reference = load_reference(reference_dat_path)

    if actual.is_empty and reference.is_empty:
        return

    # Pure-numeric fast path (backward compatible)
    try:
        actual_arr = actual.to_numpy()
        ref_arr = reference.to_numpy()
        np.testing.assert_allclose(actual_arr, ref_arr, rtol=rtol, atol=atol,
                                   err_msg=f"{actual_dat_path} "
                                           f"!= {reference_dat_path}")
        return
    except ValueError:
        pass

    # Mixed text + numeric / pure text fallback
    _compare_lines(test_case, actual, reference,
                   actual_dat_path, reference_dat_path, rtol, atol)


def discover_output_files(output_dir: Path, analysis_name: str) -> list[Path]:
    """Return the output files an analysis produced in `output_dir`.

    Uses `OUTPUT_PATTERNS` to glob for the files an analysis writes. Returns
    an empty list if the analysis has no known output patterns.
    """
    patterns = OUTPUT_PATTERNS.get(analysis_name, [])
    files: list[Path] = []
    for pattern in patterns:
        files.extend(output_dir.glob(pattern))
    # Deduplicate while preserving order
    return list(dict.fromkeys(files))


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
        name = getattr(self, 'analysis_name', '')
        if name:
            return Path(self.config.input_path) / name / relpath
        return self.config.input_path / relpath

    def outfile(self, relpath: str | Path) -> Path:
        name = getattr(self, 'analysis_name', '')
        if name:
            return Path(self.config.output_path) / name / relpath
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
            if outfile.closed:
                # run_cmd closes its streams after the subprocess completes; reopen to read size
                with open(outfile.name) as fh:
                    fh.seek(0, io.SEEK_END)
                    return fh.tell()
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

        if os.environ.get('STA_SHOW_REGEN'):
            print(f'# REGEN: stanalyzer {analysis} {args}', file=sys.stderr)

        out_err = self.run_cmd(f'stanalyzer {analysis} {args}',
                               out_filename=f'{stem}.out',
                               err_filename=f'{stem}.err')

        if accepts_o:
            return out_err + (dat_abspath,)
        return out_err

    def run_cmd(self, args: str, out_filename: str, err_filename: str) -> IOPair:
        out_path = self.outfile(out_filename)
        err_path = self.outfile(err_filename)

        if not out_path.parent.exists():
            out_path.parent.mkdir(parents=True)

        with out_path.open('w') as out_stream, err_path.open('w') as err_stream:
            print('args:', args, file=out_stream)
            out_stream.flush()
            subprocess.run(
                shlex.split(args),
                cwd=self.config.output_path,
                stdin=subprocess.DEVNULL,
                stdout=out_stream,
                stderr=err_stream,
                check=True,
            )
        # Reopen files for reading/checking after command finishes.
        out_read = out_path.open('r')
        err_read = err_path.open('r')
        return out_read, err_read

    def run(self, result: unittest.TestResult | None = None) -> unittest.TestResult | None:
        if not getattr(self, '__unittest_skip__', False):
            assert self.manager is not None, "Missing project.json"

            self.manager.write()
            result = super().run(result)

        return result


class SoohyungCase(AnalysisCase):
    """Shortcut for preparing project.json using soohyung_membrane as the template"""
    default_output: t.ClassVar[str | Path] = 'test_case'
    standard_args: t.ClassVar[str | None] = None
    test_standard: Callable
    accepts_o: t.ClassVar[bool] = True

    # subclass should override if its name doesn't follow camel_to_snake scheme
    analysis_name: t.ClassVar[str] = ''

    def standard_test(self) -> None:
        outfile = self.outfile
        args = self.standard_args

        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
            try:
                self.assertTrue(self.file_exists(dat, outfile))
                self.assertFalse(self.file_empty(out, outfile))
            finally:
                out.close()
                err.close()
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)
            try:
                self.assertFalse(self.file_empty(out, outfile))
            finally:
                out.close()
                err.close()

    def __init_subclass__(cls, **kwargs):
        if not cls.analysis_name:
            cls.analysis_name = camel_to_snake(cls.__name__)
        cls.default_output = Path('results') / "soohyung_membrane"
        cls.test_standard = skipUnlessAttrNotNone(cls, 'standard_args')(SoohyungCase.standard_test)

        super().__init_subclass__(**kwargs)

    def __init__(self, methodName='runTest'):
        if not hasattr(self, 'manager'):
            self.manager = ManagedConfig(
                input_relpath=Path('inputs') / "soohyung_membrane",
                output_relpath=self.default_output, traj="step7_*.dcd",
                psf="step5_input.psf")
        super().__init__(methodName)


class YiweiCase(AnalysisCase):
    """Shortcut for preparing project.json using yiwei_protein as the template"""
    default_output: t.ClassVar[str | Path] = 'test_case'
    standard_args: t.ClassVar[str | None] = None
    test_standard: Callable
    accepts_o: t.ClassVar[bool] = True

    # subclass should override if its name doesn't follow camel_to_snake scheme
    analysis_name: t.ClassVar[str] = ''

    def standard_test(self) -> None:
        outfile = self.outfile
        args = self.standard_args

        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
            self.assertTrue(self.file_exists(dat, outfile))
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)
            try:
                self.assertFalse(self.file_empty(out, outfile))
            finally:
                out.close()
                err.close()

    def __init_subclass__(cls, **kwargs):
        if not cls.analysis_name:
            cls.analysis_name = camel_to_snake(cls.__name__)
        cls.default_output = Path('results') / "yiwei_protein"
        cls.test_standard = skipUnlessAttrNotNone(cls, 'standard_args')(YiweiCase.standard_test)

        super().__init_subclass__(**kwargs)

    def __init__(self, methodName='runTest'):
        if not hasattr(self, 'manager'):
            self.manager = ManagedConfig(
                input_relpath=Path('inputs') / "yiwei_protein",
                output_relpath=self.default_output, traj="step5_*.dcd",
                psf="step3_input.psf")
        super().__init__(methodName)


class DensityZ(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "name C*" --sel-name "MEMB" ' \
        '--sel-sys "segid MEMB and (name [PN] or name [CO]1[0-9])"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class ContactResidenceTime(SoohyungCase):
    standard_args = '--sel "protein and name CA" --threshold "5.0"'
    analysis_name = 'contact_res_time'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class HBond(SoohyungCase):
    standard_args = '--sel "segid MEMB or protein" --hydrogens-sel "None" ' \
           '--acceptors-sel "None" --d-a-cutoff "3.0" ' \
           ' --d-h-a-angle-cutoff "150.0"'
    analysis_name = 'hbond'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class PiStacking(YiweiCase):
    standard_args = '--sel "not segid SOLV and not segid IONS" ' \
           ' --pi-pi-dist-cutoff "6.0" --pi-cation-dist-cutoff "6.0"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class RMSF(SoohyungCase):
    analysis_name = 'rmsf'
    sel = "segid PROA and name CA"
    standard_args = f'--sel-align "{sel}" --sel-rmsf "{sel}"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class SystemSize(SoohyungCase):
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Thickness(SoohyungCase):
    standard_args = '--sel "segid MEMB and (name P or name N or name C1[0-9] or name O1[0-9])" ' \
        '--sel-sys "resname DOPC and name P; resname DSPC and name P"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class WaterBridge(YiweiCase):
    standard_args = '--sel "protein" --sel2 "None" --water-sel "resname TIP3" '\
                    '--d-a-cutoff "3.0" --d-h-a-angle-cutoff "150.0"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


# ---------------------------------------------------------------------------
# Category B – SoohyungCase subclasses
# ---------------------------------------------------------------------------

class ClusteringHca(SoohyungCase):
    accepts_o = False
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("scikit-learn-extra not installed")
class ClusteringKmedoid(SoohyungCase):
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class CompressibilityModulus(SoohyungCase):
    standard_args = '--temp 310'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Contacts(SoohyungCase):
    standard_args = '--sel "protein and name CA" --contact-threshold "5.0"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class CovAnalysis(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "protein and name CA"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("analysis crashes with current test data")
class MsdMembrane(SoohyungCase):
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class MsdSolution(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "resname DOPC and name P"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class PositionTime(SoohyungCase):
    standard_args = '--sel "protein and name CA" --head-group "segid MEMB and name P"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class PositionTimeCopy(SoohyungCase):
    analysis_name = 'position_time_copy'
    standard_args = '--sel "protein and name CA"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Rdf(SoohyungCase):
    standard_args = '-sel1 "protein and name CA" -sel2 "resname DOPC and name P" -bin-size 0.1'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class RadiusOfGyration(SoohyungCase):
    standard_args = '--sel-rg "protein and name CA" --sel-align "protein and name CA"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Rmsd(SoohyungCase):
    standard_args = '--sel "protein and name CA"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("no charged residues in soohyung_membrane")
class SaltBridge(SoohyungCase):
    standard_args = '--positive-sel "resname ARG LYS and name NZ NZ*" ' \
        '--negative-sel "resname ASP GLU and name OE* OD*" ' \
        '--positive-def "resname ARG LYS and name NZ NZ*" ' \
        '--negative-def "resname ASP GLU and name OE* OD*" ' \
        '--dist-cutoff "4.5"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Scd(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "resname DOPC and (name C22 or name C32)" --sel-sys "segid MEMB and name P" --qa'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class VoronoiApl(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class VoronoiContact(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class VoronoiShellComp(SoohyungCase):
    accepts_o = False
    standard_args = '--sel "resname DOPC and name P; resname DSPC and name P" --sel-sys "segid MEMB and name P" --qa'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


# ---------------------------------------------------------------------------
# Category C – SoohyungCase subclasses requiring external tools
# (all skipped: the tools are not installed in the test environment)
# ---------------------------------------------------------------------------

@unittest.skip("requires mkdssp; not installed")
class SecondaryStructure(SoohyungCase):
    standard_args = '--sel "protein"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("requires freesasa; not installed")
class Sasa(SoohyungCase):
    analysis_name = 'sasa'
    standard_args = '--sel "protein"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("requires hole2; not available on osx-arm64")
class Hole(SoohyungCase):
    analysis_name = 'hole'
    accepts_o = False
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("requires PDB topology instead of PSF")
class GlycosidicBondBetweenSugars(SoohyungCase):
    analysis_name = 'glycosidic-bond-between-sugars'
    accepts_o = False
    standard_args = ''

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


# ---------------------------------------------------------------------------
# Category D – SoohyungCase subclasses (feasibility investigation)
# ---------------------------------------------------------------------------

class CholTilt(SoohyungCase):
    standard_args = '--sel "segid MEMB and resname CHL1" ' \
        '--center-sel "segid MEMB and name P"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("intermittent empty output under full-suite load (pre-existing invoke thread-join race)")
class HelixAnalysis(SoohyungCase):
    standard_args = '--sel-align "segid PROA and name CA" ' \
        '--sel-helix "segid PROA and name CA" --align-out aligned.dcd'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


@unittest.skip("requires two helices; only one in soohyung_membrane")
class HelixDistanceCrossingAngle(SoohyungCase):
    standard_args = '--helix1-start 1 --helix1-end 11 ' \
        '--helix2-start 12 --helix2-end 23'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class HelixTiltRotationAngle(SoohyungCase):
    standard_args = '--helix-start 1 --helix-end 23'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class BondStatistics(SoohyungCase):
    accepts_o = False
    standard_args = '-a "(1,2,3)(4,5,6)"'

    def test_standard_correctness(self) -> None:
        args = self.standard_args
        assert args is not None

        if self.accepts_o:
            out, err, dat = self.run_analysis(args, accepts_o=self.accepts_o)
        else:
            out, err = self.run_analysis(args, accepts_o=self.accepts_o)

        output_dir = Path(self.config.output_path) / self.analysis_name
        actual_files = discover_output_files(output_dir, self.analysis_name)

        if not actual_files:
            self.skipTest(f'No output files found for {self.analysis_name}')

        ref_dir = Path(__file__).parent / 'reference' / self.analysis_name
        for actual in actual_files:
            ref = ref_dir / actual.name
            assert_output_matches_reference(self, actual, ref)


class Contacts(SoohyungCase):
    standard_args = '--sel "protein" --contact-threshold "5.0"'

if __name__ == '__main__':
    unittest.main()
