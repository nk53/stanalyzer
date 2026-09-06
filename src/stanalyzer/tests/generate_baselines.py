#!/usr/bin/env python3
"""Generate golden baseline files for stanalyzer test suite.

Run standalone:
    cd src/stanalyzer/tests && python generate_baselines.py

Pre-validation: for each analysis, check that all standard_args are recognized
by the parser. Abort on any validation failure.

Generates baselines for the soohyung_membrane system (Categories A, B, C)
and the yiwei_protein system (Category Y). Use --only to regenerate a subset.
"""

import argparse
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Category definitions
# ---------------------------------------------------------------------------

CATEGORY_A = {
    'density_z': ('--sel "name C*" --sel-name "MEMB" '
                   '--sel-sys "segid MEMB and (name [PN] or name [CO]1[0-9])"'),
    'contact_res_time': '--sel "protein and name CA" --threshold "5.0"',
    'hbond': ('--sel "segid MEMB or protein" --hydrogens-sel "None" '
              '--acceptors-sel "None" --d-a-cutoff "3.0" '
              '--d-h-a-angle-cutoff "150.0"'),
    'rmsf': ('--sel-align "segid PROA and name CA" '
             '--sel-rmsf "segid PROA and name CA"'),
    'system_size': '',
    'thickness': ('--sel "segid MEMB and (name P or name N or '
                  'name C1[0-9] or name O1[0-9])" '
                  '--sel-sys "resname DOPC and name P; '
                  'resname DSPC and name P"'),
}

CATEGORY_B = {
    'salt_bridge': None,  # No charged residues in soohyung_membrane
    'contacts': '--sel "protein and name CA" --contact-threshold "5.0"',
    'voronoi_shell_comp': ('--sel "resname DOPC and name P; resname DSPC and name P" '
                           '--sel-sys "segid MEMB and name P" --qa'),
    'voronoi_contact': ('--sel "resname DOPC and name P; resname DSPC and name P" '
                        '--sel-sys "segid MEMB and name P" --qa'),
    'voronoi_apl': ('--sel "resname DOPC and name P; resname DSPC and name P" '
                    '--sel-sys "segid MEMB and name P" --qa'),
    'scd': ('--sel "resname DOPC and (name C22 or name C32)" '
            '--sel-sys "segid MEMB and name P" --qa'),
    'clustering_kmedoid': None,  # requires sklearn_extra
    'rdf': ('-sel1 "protein and name CA" -sel2 "resname DOPC and name P" '
            '-bin-size 0.1'),
    'msd_solution': '--sel "resname DOPC and name P"',
    'msd_membrane': None,  # crashes silently with exit 1, no output
    'compressibility_modulus': '--temp 310',
    'radius_of_gyration': ('--sel-rg "protein and name CA" '
                           '--sel-align "protein and name CA"'),
    'rmsd': '--sel "protein and name CA"',
    'position_time': ('--sel "protein and name CA" '
                      '--head-group "segid MEMB and name P"'),
    'position_time_copy': '--sel "protein and name CA"',
    'clustering_hca': '',
    'cov_analysis': '--sel "protein and name CA"',
}

CATEGORY_C = {
    'secondary_structure': '--sel "protein"',  # needs dssp
    'sasa': '--sel "protein"',  # needs freesasa
}

CATEGORY_Y = {
    'pi_stacking': ('--sel "not segid SOLV and not segid IONS" '
                    '--pi-pi-dist-cutoff "6.0" --pi-cation-dist-cutoff "6.0"'),
    'water_bridge': ('--sel "protein" --sel2 "None" '
                     '--water-sel "resname TIP3" --d-a-cutoff "3.0" '
                     '--d-h-a-angle-cutoff "150.0"'),
}

# Systems to process: (categories, input_dirname, traj, psf)
SYSTEMS = [
    ({**CATEGORY_A, **CATEGORY_B, **CATEGORY_C},
     'soohyung_membrane', 'step7_*.dcd', 'step5_input.psf'),
    (CATEGORY_Y, 'yiwei_protein', 'step5_*.dcd', 'step3_input.psf'),
]

# Maps analysis -> external tool dependency
TOOL_DEPS = {
    'secondary_structure': 'dssp',
    'sasa': 'freesasa',
    'clustering_kmedoid': 'sklearn_extra',
}

# ---------------------------------------------------------------------------
# Analyses that DON'T accept --out (hardcoded filenames or stdout only)
# ---------------------------------------------------------------------------

NO_OUT_ANALYSES = {
    'density_z', 'voronoi_shell_comp', 'voronoi_contact', 'voronoi_apl',
    'scd', 'msd_solution', 'msd_membrane',
}

# Output patterns per analysis — for discovering files post-run.
# Matches test_cli.py OUTPUT_PATTERNS.
OUTPUT_PATTERNS = {
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
    'density_z': ['*_nb*_*.dat', 'combined_nb*_*.dat', 'NA_*_nb*_*.dat'],
    'scd': ['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat'],
    'voronoi_shell_comp': ['ave_*_*.dat', 'time_*_*.dat', 'NA_*_*.dat'],
    'voronoi_contact': ['ave_*_*.dat', 'time_*_*.dat', 'NA_time_*_*.dat'],
    'voronoi_apl': ['ave_*_*.dat', 'time_*_*.dat'],
    'msd_solution': ['sys_com_*.dat', 'mol_com_*.dat', '*_*.dat',
                     'NA_*_*.dat', 'mol_info_*.dat'],
    'msd_membrane': ['*_sys_com_*.dat', '*_mol_com_*.dat', '*_*_*.dat',
                     'NA_*_*_*.dat', '*_mol_info_*.dat'],
    'clustering_hca': ['cluster.dat', 'cluster_representative.pdb'],
    'cov_analysis': ['corr_matrix.dat', 'eigenvalues.dat', 'eigenvectors.dat'],
    'clustering_kmedoid': ['*.dat'],
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def find_stanalyzer() -> str:
    """Find the stanalyzer CLI entry point.

    Returns a command string usable in subprocess.run(shell=True).
    Checks the same bin/ directory as the running Python first,
    then PATH, then falls back to ``python -m stanalyzer.cli.stanalyzer``.
    """
    sta = shutil.which('stanalyzer')
    if sta:
        return sta
    bin_dir = Path(sys.executable).parent
    sta_in_bin = bin_dir / 'stanalyzer'
    if sta_in_bin.exists():
        return str(sta_in_bin)
    python = sys.executable
    return f'{python} -m stanalyzer.cli.stanalyzer'


def discover_tools() -> dict[str, bool]:
    """Detect availability of optional external tools and Python modules."""
    result: dict[str, bool] = {
        'dssp': shutil.which('dssp') is not None,
        'freesasa': shutil.which('freesasa') is not None,
        'hole2': False,
    }
    try:
        import importlib
        importlib.import_module('sklearn_extra')
        result['sklearn_extra'] = True
    except ImportError:
        result['sklearn_extra'] = False
    return result


def discover_analyses() -> list[str]:
    """List analysis module names from stanalyzer.analysis directory."""
    import stanalyzer
    analysis_dir = Path(stanalyzer.__path__[0]) / 'analysis'
    analyses: list[str] = []
    for path in sorted(analysis_dir.iterdir()):
        if not path.name.endswith('.py'):
            continue
        has_main = has_parser = False
        with path.open() as f:
            for line in f:
                if line.startswith('def main'):
                    has_main = True
                elif line.startswith('def get_parser'):
                    has_parser = True
        if has_main and has_parser:
            analyses.append(path.stem)
    return analyses


def has_out_arg(analysis_name: str) -> bool:
    """Check if an analysis parser includes --out by inspecting source."""
    import stanalyzer
    source = (Path(stanalyzer.__path__[0]) / 'analysis'
              / f'{analysis_name}.py').read_text()
    return "'out'" in source and 'add_project_args' in source


def discover_output_files(output_dir: Path,
                          analysis_name: str) -> list[Path]:
    """Discover output files for an analysis in the given directory."""
    patterns = OUTPUT_PATTERNS.get(analysis_name, ['*.dat'])
    files: list[Path] = []
    for pattern in patterns:
        files.extend(output_dir.glob(pattern))
    return list(dict.fromkeys(files))


def write_project_json(output_dir: Path,
                       input_dirname: str = 'soohyung_membrane',
                       traj: str = 'step7_*.dcd',
                       psf: str = 'step5_input.psf') -> None:
    """Write project.json for the given input system.

    Writes a dict directly — matches ManagedConfig.write() pattern but
    avoids importing the Pydantic Project model for type-safety reasons.
    """
    cwd = Path.cwd().resolve()
    input_path = cwd / 'inputs' / input_dirname
    output_path = output_dir.resolve()

    project_dict: dict = {
        'title': 'Baseline Generation',
        'input_path': str(input_path),
        'output_path': str(output_path),
        'traj': traj,
        'psf': psf,
        'time_step': '1 ns',
        'scheduler': 'interactive',
    }

    out_file = output_dir / 'project.json'
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open('w') as f:
        json.dump(project_dict, f, indent=4)


# ---------------------------------------------------------------------------
# Pre-validation
# ---------------------------------------------------------------------------


def validate_analysis(analysis_name: str,
                      standard_args_str: str) -> list[str]:
    """Validate that all standard_args are recognized by the analysis parser.

    Returns a list of validation errors (empty = all good).
    Checks that each flag in standard_args exists in the parser's registered
    option strings, then does a full parse test with shlex-split args.
    """
    from importlib import import_module

    errors: list[str] = []

    try:
        module = import_module(f'stanalyzer.analysis.{analysis_name}')
    except Exception as e:
        return [f'Import failed: {e}']

    parser = module.get_parser()

    # Collect all flags the parser recognizes
    known_flags: set[str] = set()
    for action in parser._actions:
        for opt in action.option_strings:
            known_flags.add(opt)

    # Parse the standard_args_str with proper shell quoting
    if not standard_args_str.strip():
        return []  # no args to validate

    try:
        tokens = shlex.split(standard_args_str)
    except ValueError as e:
        return [f'Could not parse standard_args: {e}']

    # Walk tokens and check each flag is known
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token.startswith('-'):
            if token not in known_flags:
                errors.append(f'Unknown flag: {token}')
            # Skip the value(s) this flag consumes
            action = _find_action_for_flag(parser, token)
            if action is not None:
                nargs = action.nargs
                if nargs is None:  # consumes exactly 1
                    i += 1
                elif isinstance(nargs, int):
                    i += nargs
                elif nargs == '+':
                    # skip until next flag or end
                    i += 1
                    while i < len(tokens) and not tokens[i].startswith('-'):
                        i += 1
                    continue  # already past the last consumed token
        i += 1

    return errors


def _find_action_for_flag(parser: argparse.ArgumentParser,
                          flag: str):
    """Find the argparse action that handles a given flag string."""
    for action in parser._actions:
        if flag in action.option_strings:
            return action
    return None


# ---------------------------------------------------------------------------
# Baseline generation
# ---------------------------------------------------------------------------


def run_analysis_cmd(analysis_name: str, args: str,
                     accepts_out: bool, output_dir: Path,
                     sta_cmd: str,
                     dat_filename: str = '') -> list[Path]:
    """Run a single analysis via subprocess and return discovered output files.

    Args:
        analysis_name: Name of the analysis module.
        args: Additional CLI arguments (standard_args).
        accepts_out: Whether the analysis accepts --out.
        output_dir: The directory where analysis outputs are written.
        sta_cmd: Path to the stanalyzer CLI command.
        dat_filename: If accepts_out, the output file name to pass to --out.
    """
    out_args = ''
    if accepts_out:
        dat_name = dat_filename or f'{analysis_name}.dat'
        abs_out = str(output_dir / dat_name)
        out_args = f'--out "{abs_out}"'

    full_args = f'{sta_cmd} {analysis_name} {out_args} {args}'.strip()

    print(f'  Running: {full_args}')

    # Record existing files before run
    existing_files = set(output_dir.rglob('*'))

    result = subprocess.run(
        full_args,
        shell=True,
        capture_output=True,
        text=True,
        timeout=600,
        cwd=str(output_dir),
    )

    if result.returncode != 0:
        stderr = result.stderr.strip()
        print(f'  FAILED (exit {result.returncode}): {stderr[:200]}')
        return []

    # Discover new files (filter out directories and project.json)
    new_files = [f for f in output_dir.rglob('*')
                 if f.is_file() and f.name != 'project.json'
                 and f not in existing_files]

    return new_files


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Generate golden baseline files for stanalyzer test suite.')
    parser.add_argument('--only', nargs='+', metavar='ANALYSIS',
                        help='Only regenerate the given analyses '
                             '(across whichever system defines them).')
    cli_args = parser.parse_args()
    only_names = set(cli_args.only) if cli_args.only else None

    script_dir = Path(__file__).resolve().parent
    reference_dir = script_dir / 'reference'

    print('=' * 70)
    print('Baseline Generator for ST-Analyzer Test Suite')
    print('=' * 70)
    print(f'Script dir:  {script_dir}')
    print(f'Reference:   {reference_dir}')
    print()

    # Discover analyses and tools
    all_analyses = discover_analyses()
    tools = discover_tools()

    print(f'Found {len(all_analyses)} analysis modules')
    print(f'Tools: {tools}')
    print()

    # Build the flat list of analyses to process, carrying per-system inputs
    analyses_to_process: list[tuple[str, str, str, str, str]] = []
    # (name, args, input_dirname, traj, psf)
    for categories, input_dirname, traj, psf in SYSTEMS:
        for name, args in categories.items():
            if name not in all_analyses:
                print(f'WARNING: {name} not found in analysis modules, skipping')
                continue
            if args is None:
                print(f'SKIP {name}: incompatible with test data or missing deps')
                continue
            tool_dep = TOOL_DEPS.get(name)
            if tool_dep and not tools.get(tool_dep, False):
                print(f'SKIP {name}: requires {tool_dep} (not available)')
                continue
            if only_names is not None and name not in only_names:
                continue
            analyses_to_process.append((name, args, input_dirname, traj, psf))

    if only_names is not None:
        processed_names = {name for name, *_ in analyses_to_process}
        missing = only_names - processed_names
        if missing:
            print(f'WARNING: --only names not found in any system: '
                  f'{sorted(missing)}')

    print(f'Analyses to process: {len(analyses_to_process)}')
    print()

    sta_cmd = find_stanalyzer()
    print(f'Using stanalyzer: {sta_cmd}')
    print()

    # -----------------------------------------------------------------------
    # Phase 1: Pre-validation
    # -----------------------------------------------------------------------
    print('Phase 1: Pre-validation')
    print('-' * 70)

    validation_errors: dict[str, list[str]] = {}
    for name, args, _input_dirname, _traj, _psf in analyses_to_process:
        errors = validate_analysis(name, args)
        if errors:
            validation_errors[name] = errors
            print(f'  FAIL  {name}: {errors}')
        else:
            print(f'  OK    {name}')

    print()

    if validation_errors:
        print('PRE-VALIDATION FAILED')
        print('The following analyses have unrecognized arguments:')
        for name, errs in validation_errors.items():
            for e in errs:
                print(f'  {name}: {e}')
        print()
        print('Aborting. Fix the standard_args in this script and re-run.')
        return 1

    print('All pre-validations passed.')
    print()

    # -----------------------------------------------------------------------
    # Phase 2: Create output directory
    # -----------------------------------------------------------------------
    print('Phase 2: Setting up project configuration')
    print('-' * 70)

    # Create temp output directory (relative to cwd, which is the tests dir)
    output_relpath = Path('results') / 'baselines'
    output_dir = output_relpath.resolve()
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    print(f'  Output directory: {output_dir}')
    print()

    # -----------------------------------------------------------------------
    # Phase 3: Run analyses and generate baselines
    # -----------------------------------------------------------------------
    print('Phase 3: Generating baselines')
    print('-' * 70)

    succeeded: list[str] = []
    failed: list[str] = []
    no_output: list[str] = []

    # Group analyses by system, preserving SYSTEMS order
    systems_to_process: dict[tuple[str, str, str], list[tuple[str, str]]] = {}
    for name, args, input_dirname, traj, psf in analyses_to_process:
        key = (input_dirname, traj, psf)
        systems_to_process.setdefault(key, []).append((name, args))

    for (input_dirname, traj, psf), system_analyses in systems_to_process.items():
        print(f'\n=== System: {input_dirname} (traj={traj}, psf={psf}) ===')

        write_project_json(output_dir, input_dirname=input_dirname,
                           traj=traj, psf=psf)
        print(f'  Created project.json in {output_dir}')

        for name, args in system_analyses:
            print(f'\n--- {name} ---')

            # Clean output dir between analyses to prevent file contamination
            for item in output_dir.iterdir():
                if item.name == 'project.json':
                    continue
                if item.is_dir():
                    shutil.rmtree(item)
                else:
                    item.unlink()

            accepts_out = has_out_arg(name)

            new_files = run_analysis_cmd(name, args, accepts_out, output_dir,
                                         sta_cmd)

            if not new_files:
                print(f'  WARNING: No output files found')
                no_output.append(name)
                failed.append(name)
                continue

            # Copy to reference directory
            ref_analysis_dir = reference_dir / name
            if ref_analysis_dir.exists():
                shutil.rmtree(ref_analysis_dir)
            ref_analysis_dir.mkdir(parents=True)

            copied = 0
            for f in new_files:
                dest = ref_analysis_dir / f.name
                shutil.copy2(f, dest)
                size = dest.stat().st_size
                print(f'  Copied: {f.name} ({size} bytes)')
                copied += 1

            if copied > 0:
                succeeded.append(name)
                print(f'  SUCCESS: {copied} file(s) -> {ref_analysis_dir}')
            else:
                failed.append(name)

    # -----------------------------------------------------------------------
    # Cleanup
    # -----------------------------------------------------------------------
    print()
    print('Cleaning up temp output directory...')
    if output_dir.exists():
        shutil.rmtree(output_dir)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print()
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)
    print(f'Total processed:  {len(analyses_to_process)}')
    print(f'Succeeded:        {len(succeeded)}')
    print(f'Failed:           {len(failed)}')
    print(f'No output:        {len(no_output)}')

    if succeeded:
        print(f'\nGenerated baselines for:')
        for name in succeeded:
            ref_dir = reference_dir / name
            count = len(list(ref_dir.glob('*'))) if ref_dir.exists() else 0
            print(f'  {name}: {count} file(s)')

    if failed:
        print(f'\nFailed analyses:')
        for name in failed:
            print(f'  {name}')

    # Verify reference directory
    if ref_dirs := [d for d in reference_dir.iterdir() if d.is_dir()]:
        print(f'\nReference directories: {len(ref_dirs)}')
        for d in sorted(ref_dirs):
            files = list(d.glob('*'))
            sizes = [f.stat().st_size for f in files]
            print(f'  {d.name}/: {len(files)} file(s), '
                  f'sizes: {sizes}')
    else:
        print('\nWARNING: No reference directories created!')

    print()
    if not failed:
        print('All baselines generated successfully!')
        return 0
    else:
        print('Some baselines failed to generate.')
        return 1


if __name__ == '__main__':
    sys.exit(main())
