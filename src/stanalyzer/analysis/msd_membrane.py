#!/usr/bin/python
import argparse
import re
import typing as t
from collections.abc import Sequence

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
NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'

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

    qsplit = []
    for i in range(0, ntype):
        qsplit.append(True)

    # Process selection strings to extract information
    sel_type: list[str] = []        # selection type
    name_type: list[str] = []       # name of molecule type
    for i in range(0, ntype):
        tmps = selection[i].split()
        sel_type.append(tmps[0])   # segid/resname/moleculetype
        name_type.append(tmps[1])  # PROA/PRO*/DSPC/...

    return ProcessedArgSys(selection, ntype, qsplit, sel_type, name_type)


# Write leaflet COM
def write_leaflet_com(traj_com_sys_unwrap: list[NDFloat64],
                      framenum: int, interval: int, time_step: float,
                      nside: int, sside: list[str],
                      odir: str, suffix: str) -> None:
    for i in range(0, nside):
        sout = '#     frame       COM\n'
        sout += '#                   x           y          z\n'
        for j in range(0, framenum):
            tcom = traj_com_sys_unwrap[i][j]
            sout += f' {time_step*(interval*j+1):10.5f}'
            # sout += f' {interval*j+1:10d}'
            sout += f' {tcom[0]:10.5f} {tcom[1]:10.5f} {tcom[2]:10.5f}\n'
        # print(sout)
        sta.write_to_outfile(f'{odir}/{sside[i]}_sys_com_{suffix}.dat', sout)


# Write unwrapped COMs of individual molecule in each leaflet
def write_mol_com(traj_com_unwrap: list[NDFloat64],
                  framenum: int, interval: int, time_step: float,
                  nside: int, sside: list[str], nmol: list[int],
                  odir: str, suffix: str) -> None:
    for i in range(0, nside):
        sout = f'#  leaflet {sside[i]}\n'
        sout += '#     frame       COMs (three columns for each molecule)\n'
        sout += '#                   (x           y          z)_0 ...'
        sout += ' (x           y          z)_n: n=number of molecules\n'
        for j in range(0, framenum):
            tcom = traj_com_unwrap[i][j]
            sout += f' {time_step*(interval*j+1):10.5f}'
            sout += f' {interval*j+1:10d}'
            for k in range(0, nmol[i]):
                for m in range(0, 3):
                    sout += f' {tcom[k, m]:10.5f}'
            sout += '\n'
        # print(sout)
        sta.write_to_outfile(f'{odir}/{sside[i]}_mol_com_{suffix}.dat', sout)


# Write x,y,z-components of MSD for given molecule type in a given leaflet
def write_msd(time_step: float,
              msd: NDFloat64, taus: list[int], fout: str) -> None:
    sout = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'

    ntau = len(taus)
    for i in range(0, ntau):
        sout += f' {time_step*taus[i]:10.5f}'
        # sout += f' {taus[i]:10d}'
        for j in range(0, 3):
            sout += f' {msd[i, j]:10.5f}'
        sout += '\n'
    sta.write_to_outfile(fout, sout)


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
            sout += f' {j:10d} {name_type[jtype].strip('*')}\n'
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
                        help="Output type. outl: leaflets; outb: bilayer. Default: mda")
    parser.add_argument(
        '--qcomsys', default=False, action='store_true',
        help='If set True, write unwrapped leaflet COM time series.')
    parser.add_argument(
        '--qcommol', default=False, action='store_true',
        help='If set True, write unwrapped COMs for individual molecules.')

    return parser


def run_msd_membrane(
        sel: str, split: str, sel_sys: str, qcomsys: bool, qcommol: bool,
        psf: sta.FileRef, traj: sta.FileRefList, time_step: float | str,
        suffix: str,
        interval: int = 1,
        # lam: LeafletAssignmentMethod='mda',
        otype: OutputFileType = 'outl') -> None:
    """
    ----------
    Calculate Mean Square Displacement of COMs of selected molecule types

    MSDs are calculated separately for individual leaflets.
    Results will be obtained for leaflets or bilayer (weighted average of leaflet MSDs).
    ----------
    """

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
    print(f'MSD will be caulated every {interval} frames in lag time')

    # This is handled by ST-analyzer
    # output dir
    odir = "./"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    # bilayer recentering - should be done before any assignments
    # - center in the box (an atom)
    # - center in the box (atom group for system)
    # - unwrap to get connectd molecules
    origin = 0, 0, 0  # np.zeros([3],dtype=float) ; it did not work
    # ag_cent = u.select_atoms(sel_sys)
    ag_cent = u.atoms[[]]
    for itype in range(0, ntype_sys):
        ag_cent += u.select_atoms(selection_sys[itype])
    ag_all = u.atoms

    workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point=origin),
                transformations.center_in_box(ag_cent, point=origin),
                transformations.unwrap(ag_all)]

    u.trajectory.add_transformations(*workflow)

    # Generate ref groups for leaflet assignemnt
    if method == "zpos":
        ag_leaflet = myleaflet.assign_leaflet_zpos(u, ag_cent)
    elif method == "mda":
        ag_leaflet = myleaflet.assign_leaflet(u, ag_cent)

    print('### Generation of full atom groups for membrane molecules: START')
    # Get numbers of molecules in individual lipid types,
    # Get total number of molecules,
    # Get lipid type index for individual molecules,
    # & Generate full atom groups for leaflet COM drift correction
    #
    nmol_type_sys, nmol_sys, id_type_sys, ag_full_sys = \
        mymol.generate_full_mol_groups(
            u, ntype_sys, sel_type_sys, name_type_sys, qsplit_sys)

    print('### Generation of full atom groups for membrane molecules: DONE')

    print('### Leaflet assignment for membrane molecules: START')
    # Assign molecules to leaflet
    # id_side_sys = mymol.assign_leaflet_index(ag_full, ag_leaflet) # don't need to use
    ag_full_sys_leaflet = \
        mymol.assign_full_ag_leaflet_from_ref_leaflet(
            u, ag_full_sys, ag_leaflet)

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

    # Set arrays in use
    # For leaflets
    # smpd_arrays = [mymsd.setup_sys_mass_pos_displ_arrays(ag_leaflet[i])
    #                for i in range(nside)]
    smpd_arrays = [mymsd.setup_sys_mass_pos_displ_arrays(ag_full_sys_leaflet[i])
                   for i in range(nside)]
    suct_arrays = [mymsd.setup_unwrapped_com_traj_array(framenum)
                   for i in range(nside)]

    mass_sys       = [tmpd[0] for tmpd in smpd_arrays]  # atom masses in ind. leaflets
    tmass_sys      = [tmpd[1] for tmpd in smpd_arrays]  # mass of ind. leaflets
    # curr. atom pos. of ind. leaflets
    pos_sys        = [tmpd[2] for tmpd in smpd_arrays]
    pos_sys_prev   = [tmpd[3] for tmpd in smpd_arrays]  # prev atom pos. of ind. leaflets
    displ_sys      = [tmpd[4]    # noqa: F841
                      for tmpd in smpd_arrays]  # atom displ. of ind. leaflets
    displ_sys_com  = [tmpd[5] for tmpd in smpd_arrays]  # leaflet COM displ.

    com_sys_unwrap = [tuct[0] for tuct in suct_arrays]  # unwrappped leaflet COMs
    # traj. of unwrapped leaflet COMS
    traj_com_sys_unwrap = [tuct[1] for tuct in suct_arrays]

    # For molecules
    mpd_arrays = [mymsd.set_mass_pos_displ_arrays(nmol[i], ag[i])
                  for i in range(nside)]
    uct_arrays = [mymsd.setup_unwrapped_mol_com_traj_array(ag[i], framenum)
                  for i in range(nside)]

    mass_mol   = [tmpd[0] for tmpd in mpd_arrays]  # atom masses of ind. molecules
    tmass_mol  = [tmpd[1] for tmpd in mpd_arrays]  # mass of ind. molecules
    # current atom positions of ind. mols.
    pos        = [tmpd[2] for tmpd in mpd_arrays]
    # prev atom positions of ind. mols.
    pos_prev   = [tmpd[3] for tmpd in mpd_arrays]
    # atom displacements of ind. mols.
    displ      = [tmpd[4] for tmpd in mpd_arrays]    # noqa: F841
    # unwrapped atom positions of ind. mols.
    pos_unwrap = [tmpd[5] for tmpd in mpd_arrays]

    com_unwrap      = [tuct[0] for tuct in uct_arrays]  # unwrapped mol. COMs
    traj_com_unwrap = [tuct[1] for tuct in uct_arrays]  # traj. of unwrapped mol. COMs

    # UNWRAPPING
    print('# UNWRAP trajectories')
    # sys.exit(0)
    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]

        # do frame-wise bilayer recentering
        Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
        zref = np.zeros([2], dtype=float)
        for iside in range(0, nside):
            zref[iside] = np.mean(Lag_ref[iside].positions[:, 2])
        # translation for z-centering
        tran = 0, 0, -np.mean(zref)
        ts = transformations.translate(tran)(ts)
        ts = transformations.unwrap(ag_all)(ts)

        # get box size
        xtla, xtlb, xtlc = ts.dimensions[:3]
        box = np.array([xtla, xtlb, xtlc], dtype=float)

        # read cooridnates
        for j in range(0, nside):
            mymsd.read_coor(nmol[j], pos[j], ag[j])
            pos_sys[j] = ag_full_sys_leaflet[j].positions

            if i == 0:
                # Initialization
                # pos, unwrapped COM, and traj. of unwrapped COM for the leaflet
                mymsd.init_unwrap_sys_com(
                    pos_sys[j], mass_sys[j], tmass_sys[j], pos_sys_prev[j],
                    com_sys_unwrap[j], traj_com_sys_unwrap[j])

                # pos., unwrapped COM, and traj. of unwraped COM for ind. mol.
                mymsd.init_unwrap_mol_com(pos[j], mass_mol[j], tmass_mol[j], pos_prev[j],
                                          pos_unwrap[j], com_unwrap[j], traj_com_unwrap[j])
                # print(f'# leaflet {sside[j]}: init unwrapped com/traj done')
            else:
                # get leaflet COM displ. and update pos. (in pos_sys_prev)
                # Current pos. of leaflet is obtained inside the function
                mymsd.calculate_displ_sys_com(i, box, pos_sys[j], pos_sys_prev[j],
                                              mass_sys[j], tmass_sys[j], displ_sys_com[j])
                # print(f'# leafelt {sside[j]}: leaflet COM displ.:',displ_sys_com[j])

                # update unwrapped leaflet COM
                com_sys_unwrap[j] = com_sys_unwrap[j] + displ_sys_com[j]
                # update unwrapped leaflet COM trajectory
                np.copyto(traj_com_sys_unwrap[j][i], com_sys_unwrap[j])

                # update unwrapped mol. pos.
                mymsd.update_unwrapped_mol_pos(i, box, pos[j], pos_prev[j],
                                               pos_unwrap[j], displ_sys_com[j])

                # calculate mol. COMs & update their traj.
                mymsd.update_unwrapped_mol_com_traj(i, pos_unwrap[j], mass_mol[j],
                                                    tmass_mol[j], com_unwrap[j], traj_com_unwrap[j])

                # print(com_sys_unwrap[j])
                # print(com_unwrap[j])

    print('# UNWRAPPING TRAJ & COMS DONE')
    # sys.exit(0)

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
    taus = [interval*i for i in range(0, framenum)]
    ntau = len(taus)  # number of data points along the delay time

    # Setup msd for individual molecule types
    msd: list[NDFloat64] = []
    for i in range(0, nside):
        tmsd = mymsd.setup_msd_arrays(ntype, ntau)
        msd.append(tmsd)

    # Calculate MSD for delay times, tau in {taus}
    for i in range(0, nside):
        print(f'# leaflet {sside[i]}')
        mymsd.calculate_msd(
            taus, framenum, interval, traj_com_unwrap[i], id_type[i], msd[i])

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
