#!/usr/bin/python
import argparse
import json
import os
import re
import typing as t
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import MDAnalysis as mda
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup

import stanalyzer.cli.stanalyzer as sta
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol
from . import msd_util as mymsd

if t.TYPE_CHECKING:
    import numpy.typing as npt

ANALYSIS_NAME = 'msd_membrane'

# LeafletAssignmentMethod: t.TypeAlias = t.Literal['mda', 'zpos']
OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']
MSDEngine: t.TypeAlias = t.Literal['auto', 'direct', 'fft']
NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
COM_CACHE_VERSION = 1

# --- The following are hard set for membrane analysis
nside = 2     # up/dn
sside = ["up", "dn"]


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]


class ProcessedArgSys(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]
    sel_type: list[str]
    name_type: list[str]


class COMCacheData(t.NamedTuple):
    molecule_com: list[NDFloat64]
    system_com: list[NDFloat64] | None


def _file_identity(file_ref: sta.FileRef) -> dict[str, int | str]:
    """Return stable local-file metadata suitable for a cache identity."""
    if isinstance(file_ref, (str, os.PathLike)):
        path = Path(file_ref)
    else:
        name = getattr(file_ref, 'name', None)
        if not name:
            raise ValueError(
                "COM caching requires topology and trajectory files with paths"
            )
        path = Path(name)
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        'path': str(resolved),
        'size': stat.st_size,
        'mtime_ns': stat.st_mtime_ns,
    }


def _com_cache_metadata(
        psf: sta.FileRef,
        traj: sta.FileRefList,
        selection: list[str],
        qsplit: list[bool],
        selection_sys: list[str],
        interval: int,
        framenum: int,
        id_type: list[list[int]],
        nmol_type: list[list[int]]) -> dict[str, t.Any]:
    return {
        'version': COM_CACHE_VERSION,
        'topology': _file_identity(psf),
        'trajectories': [_file_identity(item) for item in traj],
        'selection': selection,
        'split': qsplit,
        'selection_sys': selection_sys,
        'interval': interval,
        'framenum': framenum,
        'id_type': id_type,
        'nmol_type': nmol_type,
    }


def _load_com_cache(
        cache_dir: Path,
        expected_metadata: dict[str, t.Any],
        molecule_counts: Sequence[int],
        require_system_com: bool) -> COMCacheData | None:
    metadata_path = cache_dir / 'metadata.json'
    try:
        with metadata_path.open(encoding='utf-8') as stream:
            metadata = json.load(stream)
        has_system_com = bool(metadata.pop('has_system_com', False))
        if metadata != expected_metadata:
            return None

        molecule_com = [
            np.load(cache_dir / f'{side}_mol_com.npy', mmap_mode='r')
            for side in sside
        ]
        expected_frames = expected_metadata['framenum']
        for side, values, molecule_count in zip(
                sside, molecule_com, molecule_counts, strict=True):
            if values.shape != (expected_frames, molecule_count, 3):
                raise ValueError(
                    f"cached {side} molecule COM array has shape "
                    f"{values.shape}"
                )

        system_com: list[NDFloat64] | None = None
        if require_system_com:
            if not has_system_com:
                return None
            system_com = [
                np.load(cache_dir / f'{side}_sys_com.npy', mmap_mode='r')
                for side in sside
            ]
            if any(values.shape != (expected_frames, 3)
                   for values in system_com):
                raise ValueError("cached system COM array has an invalid shape")
        return COMCacheData(molecule_com, system_com)
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


def _atomic_save_array(path: Path, values: NDFloat64) -> None:
    temporary = path.with_name(f'.{path.name}.tmp')
    with temporary.open('wb') as stream:
        np.save(stream, values, allow_pickle=False)
    os.replace(temporary, path)


def _write_com_cache(
        cache_dir: Path,
        metadata: dict[str, t.Any],
        molecule_com: Sequence[NDFloat64],
        system_com: Sequence[NDFloat64] | None) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    for side, values in zip(sside, molecule_com, strict=True):
        _atomic_save_array(cache_dir / f'{side}_mol_com.npy', values)
    if system_com is not None:
        for side, values in zip(sside, system_com, strict=True):
            _atomic_save_array(cache_dir / f'{side}_sys_com.npy', values)

    metadata_path = cache_dir / 'metadata.json'
    temporary = metadata_path.with_name(f'.{metadata_path.name}.tmp')
    stored_metadata = {
        **metadata,
        'has_system_com': system_com is not None,
    }
    with temporary.open('w', encoding='utf-8') as stream:
        json.dump(stored_metadata, stream, indent=2, sort_keys=True)
        stream.write('\n')
    os.replace(temporary, metadata_path)


def assign_leaflet_zpos_fast(atomgroup):
    """Split an existing atom group by its z midpoint without global selects."""
    positions = atomgroup.positions
    z_positions = positions[:, 2]
    z_center = 0.5 * (
        np.min(z_positions)
        + np.max(z_positions)
    )
    return [
        atomgroup[z_positions > z_center],
        atomgroup[z_positions < z_center],
        atomgroup[z_positions == z_center],
    ]


def process_args(sel: str, split_to_mol: str | None) -> ProcessedArgs:
    """
    ----------
    Process arguments
    ----------
    """

    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()
    if any(not item for item in selection):
        raise ValueError("sel contains an empty molecule selection")

    # in case split_to_mol is None
    if not split_to_mol:
        nsplit = 0
    # normal process of split_to_mol
    else:
        split = re.split(';|,', f'{split_to_mol:s}')
        nsplit = len(split)
        for i in range(0, nsplit):
            split[i] = split[i].strip()

    qsplit = []
    if nsplit < ntype:  # add more qsplit options
        for i in range(0, nsplit):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)
        for i in range(nsplit, ntype):
            qsplit.append(True)  # default values
    else:  # get split up to ntype
        qsplit = []
        for i in range(0, ntype):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)

    return ProcessedArgs(selection, ntype, qsplit)


def process_arg_sys(sel: str) -> ProcessedArgSys:
    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()
    if any(not item for item in selection):
        raise ValueError("sel_sys contains an empty membrane selection")

    qsplit = []
    for i in range(0, ntype):
        qsplit.append(True)

    # Process selection strings to extract information
    sel_type: list[str] = []        # selection type
    name_type: list[str] = []       # name of molecule type
    for i in range(0, ntype):
        tmps = selection[i].split()
        if len(tmps) < 2:
            raise ValueError(
                "Each sel_sys entry must begin with a selection type "
                "and molecule name, for example 'resname DMPC'"
            )
        sel_type.append(tmps[0])   # segid/resname/moleculetype
        name_type.append(tmps[1])  # PROA/PRO*/DSPC/...

    return ProcessedArgSys(selection, ntype, qsplit, sel_type, name_type)


# Write leaflet COM
def write_leaflet_com(traj_com_sys_unwrap: list[NDFloat64],
                      framenum: int, interval: int, time_step: float,
                      nside: int, sside: list[str],
                      odir: str, suffix: str) -> None:
    for i in range(0, nside):
        fout = f'{odir}/{sside[i]}_sys_com_{suffix}.dat'
        header = '#     frame       COM\n'
        header += '#                   x           y          z\n'
        stream = sta.write_to_outfile(fout, header, close=False)
        try:
            for j in range(0, framenum):
                tcom = traj_com_sys_unwrap[i][j]
                stream.write(
                    f' {time_step*(interval*j+1):10.5f}'
                    f' {tcom[0]:10.5f} {tcom[1]:10.5f}'
                    f' {tcom[2]:10.5f}\n'
                )
        finally:
            stream.close()


# Write unwrapped COMs of individual molecule in each leaflet
def write_mol_com(traj_com_unwrap: list[NDFloat64],
                  framenum: int, interval: int, time_step: float,
                  nside: int, sside: list[str], nmol: list[int],
                  odir: str, suffix: str) -> None:
    for i in range(0, nside):
        fout = f'{odir}/{sside[i]}_mol_com_{suffix}.dat'
        header = f'#  leaflet {sside[i]}\n'
        header += '#     frame       COMs (three columns for each molecule)\n'
        header += '#                   (x           y          z)_0 ...'
        header += ' (x           y          z)_n: n=number of molecules\n'
        stream = sta.write_to_outfile(fout, header, close=False)
        try:
            for j in range(0, framenum):
                tcom = traj_com_unwrap[i][j]
                row = [
                    f' {time_step*(interval*j+1):10.5f}',
                    f' {interval*j+1:10d}',
                ]
                row.extend(
                    f' {value:10.5f}'
                    for value in tcom[:nmol[i]].flat
                )
                stream.write(''.join(row))
                stream.write('\n')
        finally:
            stream.close()


# Write x,y,z-components of MSD for given molecule type in a given leaflet
def write_msd(time_step: float,
              msd: NDFloat64, taus: list[int], fout: str) -> None:
    header = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'
    stream = sta.write_to_outfile(fout, header, close=False)
    try:
        for tau, values in zip(taus, msd, strict=True):
            stream.write(
                f' {time_step*tau:10.5f}'
                f' {values[0]:10.5f}'
                f' {values[1]:10.5f}'
                f' {values[2]:10.5f}\n'
            )
    finally:
        stream.close()


# Write MSD outputs for leaflets
def write_msd_outputs_leaflet(time_step: float,
                              msd: Sequence[NDFloat64], taus: list[int],
                              nside: int, sside: list[str],
                              ntype: int, name_type: list[str], numb_type: list[list[int]],
                              odir: str, suffix: str) -> None:
    for i in range(0, nside):
        side = sside[i]
        tnumb_type = numb_type[i]

        for j in range(0, ntype):
            tnamej = name_type[j].strip('*')
            print(f'# Write MSDs for {tnamej} in leafelt, {side}')

            if tnumb_type[j] == 0:  # NA
                fout = f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat'
                sout = f'no {tnamej} in {side} leaflet.\n'
                sta.write_to_outfile(fout, sout)
            else:
                fout = f'{odir}/{side}_{tnamej.lower()}_{suffix}.dat'
                write_msd(time_step,
                          msd[i][j], taus, fout)


# Write MSD outputs for bilayer
def write_msd_outputs_bilayer(time_step: float,
                              bmsd: NDFloat64, taus: list[int],
                              ntype: int, name_type: list[str], numb_type: list[int],
                              odir: str, suffix: str) -> None:
    side = 'bilayer'
    for j in range(0, ntype):
        tnamej = name_type[j].strip('*')
        print(f'# Write MSDs for {tnamej} in bilayer')

        if numb_type[j] == 0:  # NA
            fout = f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat'
            sout = f'no {tnamej} in {side}.\n'
            sta.write_to_outfile(fout, sout)
        else:
            fout = f'{odir}/{side}_{tnamej.lower()}_{suffix}.dat'
            write_msd(time_step,
                      bmsd[j], taus, fout)


def write_mol_info(nside: int, sside: list[str], name_type: list[str],
                   nmol: list[int], id_type: list[list[int]],
                   odir: str, suffix: str) -> None:
    # write molecule info: molecule number name
    for i in range(0, nside):
        side = sside[i]
        sout = '# mol.index     name\n'
        # loop over molecules
        for j in range(0, nmol[i]):
            jtype = id_type[i][j]  # molecule type index of molecule, j
            sout += f" {j:10d} {name_type[jtype].strip('*')}\n"
        fout = f'{odir}/{side}_mol_info_{suffix}.dat'
        sta.write_to_outfile(fout, sout)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval', 'time_step')
    parser.add_argument(
        '--sel', metavar='selection', default='',
        help='Selection for individual molecule. MSD will be calculated '
             'for COMs of individual molecules.')
    parser.add_argument(
        '--split', default='',
        help='Y/N. If N, the atom group for the selection is considered as '
             'that for a single molecule. If Y, the atom group is further split '
             'to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument(
        '--sel-sys', metavar='selection', default='',
        help='Selection for system atom groups for leaflet COM drift '
             'correction and bilayer recentering.')
    # ## Commented out: Activate when non-planar bilayers are supported.
    # parser.add_argument('--lam', metavar='OPT', default='mda', choices=['mda', 'zpos'],
    #                     help='Leleat assignment method. mda: MDAnalysis.analysis.leaflet; '
    #                     'zpos: Z-position. Default: mda')
    parser.add_argument('--suffix', type=str, default='0',
                        help='Suffix  to output file(s)')
    parser.add_argument('--otype', metavar='OPT', default='outl', choices=['outl', 'outb'],
                        help="Output type. outl: leaflets; outb: bilayer. Default: outl")
    parser.add_argument(
        '--qcomsys', default=False, action='store_true',
        help='If set True, write unwrapped leaflet COM time series.')
    parser.add_argument(
        '--qcommol', default=False, action='store_true',
        help='If set True, write unwrapped COMs for individual molecules.')
    parser.add_argument(
        '--max-lag-frames', type=int, default=None,
        help='Maximum lag in analyzed frames. Default: all available lags.')
    parser.add_argument(
        '--msd-engine', choices=['auto', 'direct', 'fft'], default='auto',
        help='MSD algorithm. auto selects from frame and lag counts; direct '
             'preserves legacy reduction order; fft is faster for sufficiently '
             'long trajectories. Default: auto.')
    parser.add_argument(
        '--fft-chunk-molecules', type=int, default=256,
        help='Maximum molecules transformed together by the FFT engine. '
             'Smaller values reduce peak memory. Default: 256.')
    parser.add_argument(
        '--fft-workers', type=int, default=1,
        help='Worker threads used by SciPy FFT transforms. Default: 1.')
    parser.add_argument(
        '--com-cache', type=str, default=None,
        help='Directory used to reuse memory-mapped molecular COM trajectories.')
    parser.add_argument(
        '--refresh-com-cache', default=False, action='store_true',
        help='Ignore and replace an existing COM cache.')

    return parser


def run_msd_membrane(
        sel: str, split: str, sel_sys: str, qcomsys: bool, qcommol: bool,
        psf: sta.FileRef, traj: sta.FileRefList, time_step: float | str,
        suffix: str, interval: int = 1,
        # lam: LeafletAssignmentMethod='mda',
        otype: OutputFileType = 'outl',
        max_lag_frames: int | None = None,
        msd_engine: MSDEngine = 'auto',
        fft_chunk_molecules: int = 256,
        fft_workers: int = 1,
        com_cache: str | None = None,
        refresh_com_cache: bool = False) -> None:
    """
    ----------
    Calculate Mean Square Displacement of COMs of selected molecule types

    MSDs are calculated separately for individual leaflets.
    Results will be obtained for leaflets or bilayer (weighted average of leaflet MSDs).
    ----------
    """
    if not sel.strip():
        raise ValueError("sel must contain at least one molecule selection")
    if not sel_sys.strip():
        raise ValueError(
            "sel_sys must contain at least one membrane reference selection"
        )
    if interval < 1:
        raise ValueError("interval must be at least 1")
    if max_lag_frames is not None and max_lag_frames < 0:
        raise ValueError("max_lag_frames must be nonnegative")
    if fft_chunk_molecules < 1:
        raise ValueError("fft_chunk_molecules must be positive")
    if fft_workers < 1:
        raise ValueError("fft_workers must be positive")

    # process arguments
    selection, ntype, qsplit = process_args(sel, split)
    # process sys arguments - To generate full atom groups for individual molecules
    selection_sys, ntype_sys, qsplit_sys, sel_type_sys, name_type_sys \
        = process_arg_sys(sel_sys)

    # method=lam
    method: t.Literal['zpos', 'mda'] = 'zpos'  # use this option for plana bilayers
    outtype = otype
    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])
    if not np.isfinite(time_step) or time_step <= 0:
        raise ValueError("time_step must be finite and positive")

    # print summary of arguments
    for i in range(0, ntype):
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    print('Suffix to output files', suffix)
    print(f'Bilayer is recentered at z = 0 using {sel_sys}')
    if method == 'zpos':
        print('Leaflets are assigned based on z-positions')
    elif method == 'mda':
        print('Leaflets are assgined using a modified version of MDA LeafletFinder')
    print('Output type', outtype)
    print('Writing unwrapped COM of leaflets:', qcomsys)
    print('Writing unwrapped COM of individual molecules:', qcommol)
    print(f'MSD will be calculated every {interval} frames in lag time')
    print('Requested MSD calculation engine:', msd_engine)
    print('FFT worker threads:', fft_workers)
    print('COM cache:', com_cache if com_cache else 'disabled')

    # This is handled by ST-analyzer
    # output dir
    odir = "./"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    n_trajectory_frames = u.trajectory.n_frames
    if n_trajectory_frames < 1:
        raise ValueError("trajectory contains no frames")
    # Preserve the original sampling contract: analyze exactly
    # floor(n_frames / interval) frames and drop a trailing remainder.
    framenum = int(n_trajectory_frames / interval)
    frame_indices = [
        interval * frame
        for frame in range(framenum)
    ]
    if framenum < 1:
        raise ValueError(
            "interval is larger than the number of trajectory frames"
        )

    # bilayer recentering - should be done before any assignments
    # - center in the box (an atom)
    # - center in the box (atom group for system)
    # - unwrap to get connectd molecules
    origin = 0, 0, 0  # np.zeros([3],dtype=float) ; it did not work
    # ag_cent = u.select_atoms(sel_sys)
    ag_cent = u.atoms[[]]
    for itype in range(0, ntype_sys):
        ag_cent += u.select_atoms(selection_sys[itype])
    if len(ag_cent) == 0:
        raise ValueError(
            f"No atoms found for sel_sys: {sel_sys}"
        )

    print('### Generation of full atom groups for membrane molecules: START')
    # Generate complete membrane molecules before installing the trajectory
    # transformations. The centering selection may contain only reference
    # atoms, but unwrapping requires complete bonded molecules.
    nmol_type_sys, nmol_sys, id_type_sys, ag_full_sys = \
        mymol.generate_full_mol_groups(
            u, ntype_sys, sel_type_sys, name_type_sys, qsplit_sys)

    ag_membrane = u.atoms[[]]
    for molecule in ag_full_sys:
        ag_membrane += molecule

    print('### Generation of full atom groups for membrane molecules: DONE')

    try:
        bonded_topology = len(u.bonds) > 0
    except mda.exceptions.NoDataError:
        bonded_topology = False

    if not bonded_topology:
        # Coordinate-only topologies such as PDB do not necessarily contain
        # bonds. Guess bonds only for the membrane so MDAnalysis can unwrap
        # its molecules without imposing guessed connectivity on the rest of
        # the system.
        print(
            "Topology has no bonds; guessing membrane bonds for unwrapping"
        )
        ag_membrane.guess_bonds()
    if bonded_topology:
        workflow = [
            transformations.center_in_box(
                AtomGroup([ag_cent[0]]),
                point=origin,
            ),
            transformations.center_in_box(ag_cent, point=origin),
            # Only membrane coordinates contribute to this analysis. Limiting
            # unwrapping to them avoids fragment and box work for solvent,
            # ions, and other unrelated atoms.
            transformations.unwrap(ag_membrane),
        ]
        u.trajectory.add_transformations(*workflow)
    else:
        origin_array = np.asarray(origin, dtype=float)
        unwrap_membrane = transformations.unwrap(ag_membrane)

        def prepare_membrane_frame(ts):
            """Center and unwrap a membrane from a bondless topology."""
            membrane_center = ag_cent.center_of_geometry()
            ts.positions += origin_array - membrane_center
            return unwrap_membrane(ts)

        u.trajectory.add_transformations(prepare_membrane_frame)

    # Generate ref groups for leaflet assignemnt
    if method == "zpos":
        if bonded_topology:
            ag_leaflet = myleaflet.assign_leaflet_zpos(u, ag_cent)
        else:
            ag_leaflet = assign_leaflet_zpos_fast(ag_cent)
    elif method == "mda":
        ag_leaflet = myleaflet.assign_leaflet(u, ag_cent)

    print('### Leaflet assignment for membrane molecules: START')
    # Assign molecules to leaflet
    # id_side_sys = mymol.assign_leaflet_index(ag_full, ag_leaflet) # don't need to use
    ag_full_sys_leaflet = \
        mymol.assign_full_ag_leaflet_from_ref_leaflet(
            u, ag_full_sys, ag_leaflet[:nside])

    print('### Leaflet assignemnt for membrane molecules: DONE')
    # sys.exit(0)

    # LEAFLETs are assigned in this stage and will not be altered.
    print('### Leaflet assignment for molecule types subject to MSD calculation: START')
    # Generate molecule groups in leaflets
    name_type, nmol_type, nmol, id_type, ag =\
        mymol.generate_mol_groups_memb(
            u, nside, ntype, selection, qsplit, sside, method)
    print('### Leaflet assignment for molecule types subject to MSD calculation: DONE')

    # get bilayer nmol_type
    nmol_type0 = np.sum(nmol_type, axis=0)

    cache_dir = Path(com_cache).expanduser() if com_cache else None
    cache_metadata: dict[str, t.Any] | None = None
    cache_data: COMCacheData | None = None
    if cache_dir is not None:
        cache_metadata = _com_cache_metadata(
            psf, traj, selection, qsplit, selection_sys,
            interval, framenum, id_type, nmol_type,
        )
        if not refresh_com_cache:
            cache_data = _load_com_cache(
                cache_dir,
                cache_metadata,
                nmol,
                require_system_com=qcomsys,
            )
    if cache_data is None:
        print('# COM cache miss; preprocessing trajectory')
    else:
        print('# COM cache hit; reusing molecular COM trajectories')

    # Set arrays in use
    # For leaflets
    # smpd_arrays = [mymsd.setup_sys_mass_pos_displ_arrays(ag_leaflet[i])
    #                for i in range(nside)]
    smpd_arrays = [mymsd.setup_sys_mass_pos_displ_arrays(ag_full_sys_leaflet[i])
                   for i in range(nside)]
    sys_com_frame_count = framenum if qcomsys else 1
    suct_arrays = [mymsd.setup_unwrapped_com_traj_array(sys_com_frame_count)
                   for i in range(nside)]

    mass_sys       = [tmpd[0] for tmpd in smpd_arrays]  # atom masses in ind. leaflets
    tmass_sys      = [tmpd[1] for tmpd in smpd_arrays]  # mass of ind. leaflets
    # curr. atom pos. of ind. leaflets
    pos_sys        = [tmpd[2] for tmpd in smpd_arrays]
    pos_sys_prev   = [tmpd[3] for tmpd in smpd_arrays]  # prev atom pos. of ind. leaflets
    displ_sys      = [tmpd[4]    # noqa: F841
                      for tmpd in smpd_arrays]  # atom displ. of ind. leaflets
    displ_sys_com  = [tmpd[5] for tmpd in smpd_arrays]  # leaflet COM displ.
    sys_atom_indices = [
        leaflet.indices
        for leaflet in ag_full_sys_leaflet
    ]
    sys_coordinate_buffers = [
        np.empty((len(indices), 3), dtype=u.trajectory.ts.positions.dtype)
        for indices in sys_atom_indices
    ]

    com_sys_unwrap = [tuct[0] for tuct in suct_arrays]  # unwrappped leaflet COMs
    # traj. of unwrapped leaflet COMS
    traj_com_sys_unwrap = [tuct[1] for tuct in suct_arrays]

    # Precompute molecule indices, masses, boundaries, and contiguous working
    # arrays once. Frame processing then gathers coordinates once per leaflet.
    packed_molecules = [
        mymsd.setup_packed_molecule_data(
            ag[i],
            framenum,
        )
        for i in range(nside)
    ]
    traj_com_unwrap = [
        data.traj_com_unwrap
        for data in packed_molecules
    ]

    # UNWRAPPING
    print('# UNWRAP trajectories')
    # sys.exit(0)
    progress_stride = max(1, framenum // 100)
    frame_iterator: t.Iterable[tuple[int, t.Any]]
    if cache_data is not None:
        frame_iterator = iter(())
    elif bonded_topology:
        frame_iterator = (
            (i, u.trajectory[interval * i])
            for i in range(framenum)
        )
    else:
        frame_iterator = enumerate(u.trajectory[::interval])
    for i, ts in frame_iterator:
        frame_index = frame_indices[i]
        #  ct=(cnt-1)+dt*(i+1) # in ns
        if (
            i == 0
            or i == framenum - 1
            or i % progress_stride == 0
        ):
            print(
                f'# processing frame {frame_index + 1}/'
                f'{n_trajectory_frames}'
            )
        # do frame-wise bilayer recentering
        if bonded_topology:
            Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
        else:
            Lag_ref = assign_leaflet_zpos_fast(ag_cent)
        zref = np.zeros([2], dtype=float)
        for iside in range(0, nside):
            zref[iside] = np.mean(Lag_ref[iside].positions[:, 2])
        # The installed trajectory transformation has already unwrapped the
        # membrane. Apply the frame-specific z recentering directly so the
        # complete transformation stack is not constructed and run twice.
        ts.positions[:, 2] -= np.mean(zref)

        # Retain all six unit-cell values so triclinic boxes can use the
        # correct minimum-image convention.
        box = np.asarray(ts.dimensions, dtype=float)
        if (
            box.shape != (6,)
            or not np.isfinite(box).all()
            or np.any(box[:3] <= 0)
        ):
            raise ValueError(
                f"Invalid periodic box dimensions at frame {frame_index}"
            )

        # read cooridnates
        for j in range(0, nside):
            mymsd.read_packed_coordinates(
                ts.positions,
                packed_molecules[j],
            )
            np.take(
                ts.positions,
                sys_atom_indices[j],
                axis=0,
                out=sys_coordinate_buffers[j],
            )
            np.copyto(pos_sys[j], sys_coordinate_buffers[j])

            if i == 0:
                # Initialization
                # pos, unwrapped COM, and traj. of unwrapped COM for the leaflet
                mymsd.init_unwrap_sys_com(
                    pos_sys[j], mass_sys[j], tmass_sys[j], pos_sys_prev[j],
                    com_sys_unwrap[j], traj_com_sys_unwrap[j])

                # pos., unwrapped COM, and traj. of unwraped COM for ind. mol.
                mymsd.init_packed_molecule_com(
                    packed_molecules[j]
                )
                # print(f'# leaflet {sside[j]}: init unwrapped com/traj done')
            else:
                # get leaflet COM displ. and update pos. (in pos_sys_prev)
                # Current pos. of leaflet is obtained inside the function
                mymsd.calculate_displ_sys_com(
                    i, box, pos_sys[j], pos_sys_prev[j],
                    mass_sys[j], tmass_sys[j], displ_sys_com[j],
                    scratch=displ_sys[j],
                )
                # print(f'# leafelt {sside[j]}: leaflet COM displ.:',displ_sys_com[j])

                # update unwrapped leaflet COM
                com_sys_unwrap[j] = com_sys_unwrap[j] + displ_sys_com[j]
                if qcomsys:
                    # update the optional unwrapped leaflet COM trajectory
                    np.copyto(traj_com_sys_unwrap[j][i], com_sys_unwrap[j])

                # update unwrapped mol. pos.
                mymsd.update_packed_molecule_com(
                    i,
                    box,
                    packed_molecules[j],
                    displ_sys_com[j],
                )

                # print(com_sys_unwrap[j])
                # print(com_unwrap[j])

    print('# UNWRAPPING TRAJ & COMS DONE')
    # sys.exit(0)

    if cache_data is not None:
        traj_com_unwrap = cache_data.molecule_com
        if cache_data.system_com is not None:
            traj_com_sys_unwrap = cache_data.system_com
        del packed_molecules
    elif cache_dir is not None and cache_metadata is not None:
        _write_com_cache(
            cache_dir,
            cache_metadata,
            traj_com_unwrap,
            traj_com_sys_unwrap if qcomsys else None,
        )
        print(f'# Wrote COM cache: {cache_dir}')

    if qcomsys:
        print('# Write unwrapped leaflet COMs')
        write_leaflet_com(traj_com_sys_unwrap, framenum, interval, time_step,
                          nside, sside, odir, suffix)

    if qcommol:
        print('# Write unwrapped mol. COMs in each leaflet')
        write_mol_com(traj_com_unwrap, framenum, interval, time_step,
                      nside, sside, nmol, odir, suffix)
        print('# Write mol. info for each leaflet: (mol. number,type)')
        write_mol_info(nside, sside, name_type, nmol, id_type, odir, suffix)

    print('# MSD calculations')
    # Loop over delay times with given interval
    maximum_lag = framenum - 1
    if max_lag_frames is not None:
        maximum_lag = min(maximum_lag, max_lag_frames)
    taus = [interval*i for i in range(0, maximum_lag + 1)]
    ntau = len(taus)  # number of data points along the delay time
    selected_engine: t.Literal['direct', 'fft']
    if msd_engine == 'auto':
        selected_engine = mymsd.select_msd_engine(
            framenum,
            range(maximum_lag + 1),
        )
    else:
        selected_engine = msd_engine
    print('Selected MSD calculation engine:', selected_engine)

    # Setup msd for individual molecule types
    msd: list[NDFloat64] = []
    for i in range(0, nside):
        tmsd = mymsd.setup_msd_arrays(ntype, ntau)
        msd.append(tmsd)

    # Calculate MSD for delay times, tau in {taus}
    for i in range(0, nside):
        print(f'# leaflet {sside[i]}')
        if selected_engine == 'direct':
            mymsd.calculate_msd(
                taus, framenum, interval,
                traj_com_unwrap[i], id_type[i], msd[i])
        elif selected_engine == 'fft':
            mymsd.calculate_msd_fft(
                taus, framenum, interval,
                traj_com_unwrap[i], id_type[i], msd[i],
                chunk_size=fft_chunk_molecules,
                workers=fft_workers)
        else:
            raise ValueError(f"Unsupported MSD engine: {selected_engine}")

        for molecule_type in range(ntype):
            if nmol_type[i][molecule_type] == 0:
                continue
            values = msd[i][molecule_type]
            if not np.isfinite(values).all():
                raise RuntimeError(
                    f"Non-finite MSD values for {name_type[molecule_type]} "
                    f"in the {sside[i]} leaflet"
                )
            if np.any(values < 0):
                raise RuntimeError(
                    f"Negative MSD values for {name_type[molecule_type]} "
                    f"in the {sside[i]} leaflet"
                )
            if not np.array_equal(
                values[0],
                np.zeros(3, dtype=values.dtype),
            ):
                raise RuntimeError(
                    f"MSD at zero lag is not zero for "
                    f"{name_type[molecule_type]} in the {sside[i]} leaflet"
                )

    # Write MSD outputs
    if outtype == 'outb':
        # calculate MSD over bilayers
        bmsd = mymsd.calculate_msd_bilayer(msd, nside, ntype, ntau, nmol_type)
        write_msd_outputs_bilayer(time_step,
                                  bmsd, taus, ntype, name_type, nmol_type0, odir, suffix)
    elif outtype == 'outl':
        write_msd_outputs_leaflet(time_step,
                                  msd, taus, nside, sside,
                                  ntype, name_type, nmol_type, odir, suffix)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_msd_membrane(**settings)


if __name__ == '__main__':
    main()
