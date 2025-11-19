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

NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'

# --- The following are hard set for membrane analysis
nside = 2     # up/dn
sside = ["up", "dn"]


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]
    sel_sys: str


def process_args(sel: str, split_str: str, sel_sys: str) -> ProcessedArgs:
    """
    ----------
    Process arguments
    ----------
    """

    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()
    split = re.split(';|,', f'{split_str:s}')
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

    sel_sys = f'{sel_sys:s}'
    return ProcessedArgs(selection, ntype, qsplit, sel_sys)


# Write system COM
def write_com_sys(traj_com_sys_unwrap: list[NDFloat64], nside: int,
                  framenum: int, interval: int, odir: str, sside: list[str]) -> None:
    for i in range(0, nside):
        sout = '#     frame       COM\n'
        sout += '#                   x           y          z\n'
        for j in range(0, framenum):
            tcom = traj_com_sys_unwrap[i][j]
            sout += f' {interval*j:10d} {tcom[0]:10.5f} {tcom[1]:10.5f} {tcom[2]:10.5f}\n'
        # print(sout)
        sta.write_to_outfile(f'{odir}/{sside[i]}_com_sys_unwrapped.dat', sout)


# Write unwrapped COMs of individual molecule
def write_com_mol(traj_com_unwrap: list[NDFloat64], nside: int, nmol:
                  list[int], framenum: int, interval: int, odir: str,
                  sside: list[str]) -> None:
    for i in range(0, nside):
        sout = f'#  leaflet {sside[i]}\n'
        sout += '#     frame       COM\n'
        sout += '#                   x           y          z\n'
        for j in range(0, framenum):
            tcom = traj_com_unwrap[i][j]
            sout += f' {interval*j:10d}'
            for k in range(0, nmol[i]):
                for m in range(0, 3):
                    sout += f' {tcom[k, m]:10.5f}'
            sout += '\n'
        # print(sout)
        sta.write_to_outfile(f'{odir}/{sside[i]}_com_mol_unwrapped.dat', sout)


# Write x,y,z-components of MSD for given molecule type in a given leaflet
def write_msd(name: str, msd: NDFloat64, taus: list[int], odir: str, side: str) -> None:
    sout = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'

    ntau = len(taus)
    for i in range(0, ntau):
        sout += f' {taus[i]:10.5f}'
        for j in range(0, 3):
            sout += f' {msd[i, j]:10.5f}'
        sout += '\n'

    sta.write_to_outfile(f'{odir}/{side}_{name.lower()}_msd.dat', sout)


# Write MSD outputs for leaflets
def write_msd_outputs_leaflet(
        msd: Sequence[NDFloat64], taus: list[int], nside: int, ntype: int,
        name_type: list[str], odir: str, sside: list[str]) -> None:
    for i in range(0, nside):
        side = sside[i]

        for j in range(0, ntype):
            name = name_type[j]

            print(f'# Write MSDs for {name} in leafelt, {side}')
            write_msd(name, msd[i][j], taus, odir, side)


# Write MSD outputs for bilayer
def write_msd_outputs_bilayer(bmsd: NDFloat64, taus: list[int], ntype: int,
                              name_type: list[str], odir: str) -> None:
    side = "bilayer"

    for i in range(0, ntype):
        name = name_type[i]

        print(f'# Write MSDs for {name} in bilayer')
        write_msd(name, bmsd[i], taus, odir, side)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval')
    parser.add_argument(
        '--center', default=False, action='store_true',
        help='Perform membrane centering.')
    parser.add_argument(
        '--sel', metavar='selection',
        help='Selection for individual molecule. MSD will be calculated '
             'for COMs of individual molecules.')
    parser.add_argument(
        '--split',
        help='Y/N. If N, the atom group for the selection is considered as '
             'that for a single molecule. If Y, the atom group is further split '
             'to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument(
        '--sel-sys', metavar='selection',
        help='Selection for system atom groups for leaflet COM drift '
             'correction and bilayer recentering.')
    parser.add_argument(
        '--qz', action='store_true', default=False,
        help='Z-position based leaflet assigment: Default is false. Maybe '
             'useful when selections are minor components of the bilayer')
    parser.add_argument(
        '--qb', action='store_true', default=False,
        help='Do bilayer analysis if provided: only for sym. bilayers')
    parser.add_argument(
        '--qcomsys', default=False, action='store_true',
        help='When provided, write unwrapped COM trajectories for leaflets.')
    parser.add_argument(
        '--qcommol', default=False, action='store_true',
        help='When provided, write unwrapped COM trajectoreis for molecules.')

    return parser


def run_msd_membrane(
        sel: str, split: str, sel_sys: str, qz: bool, qb: bool, qcomsys: bool, qcommol: bool,
        psf: sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
    """
    ----------
    Calculate Mean Square Displacement of COMs of selected molecule types

    MSDs are calculated separately for individual leaflets.
    Results will be obtained for leaflets or bilayer (weighted average of leaflet MSDs).
    ----------
    """

    # process arguments
    selection, ntype, qsplit, sel_sys = process_args(sel, split, sel_sys)

    # print summary of arguments
    for i in range(0, ntype):
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    if qz:
        print('Leaflets are assigned based on z-positions')
    else:
        print('Leaflets are assgined using a modified version of MDA LeafletFinder')
    print('Writing unwrapped COM of leaflets:', qcomsys)
    print('Writing unwrapped COM of individual molecules:', qcommol)
    print(f'MSD will be caulated every {interval} frames in delay time')
    print(f'Bilayer is recentered at z = 0 using {sel_sys}:', center)

    odir = "msd"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    # if center: bilayer recentering - should be done before any assignments
    # - center in the box (an atom)
    # - center in the box (atom group for system)
    # - unwrap to get connectd molecules
    # Taken from
    if center:
        origin = 0, 0, 0  # np.zeros([3],dtype=float) ; it did not work
        ag_cent = u.select_atoms(sel_sys)
        ag_all = u.atoms

        workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point=origin),
                    transformations.center_in_box(ag_cent, point=origin),
                    transformations.unwrap(ag_all)]

        u.trajectory.add_transformations(*workflow)

    # LEAFLETs are assigned in this stage and will not be altered.
    # Generate molecule groups
    name_type, nmol_type, nmol, id_type, ag =\
        mymol.generate_mol_groups_memb(
            u, nside, ntype, selection, qsplit, sside, qz)

    # Generate system groups: required for leaflet COM drift correction
    ag_sys = mymol.generate_sys_groups(u, sel_sys, qz)

    # Set arrays in use
    mpd_arrays = [mymsd.set_mass_pos_displ_arrays(nmol[i], ag[i])
                  for i in range(nside)]
    ctc_unwrap = [mymsd.setup_unwrapped_mol_com_traj_array(ag[i], framenum)
                  for i in range(nside)]
    smpd_arrays = [mymsd.setup_sys_mass_pos_displ_arrays(ag_sys[i])
                   for i in range(nside)]
    ctc_sys_unwrap = [mymsd.setup_unwrapped_com_traj_array(framenum)
                      for _ in range(nside)]

    # atom masses and total mass of ind. molecules
    mass_mol   = [side[0] for side in mpd_arrays]
    tmass_mol  = [side[1] for side in mpd_arrays]

    # current and previous atom positions in ind. molecules
    pos        = [side[2] for side in mpd_arrays]
    pos_prev   = [side[3] for side in mpd_arrays]

    # displacements and unwrapped position of atoms
    displ      = [side[4] for side in mpd_arrays]
    pos_unwrap = [side[5] for side in mpd_arrays]

    # unwrapped COM and its traj of ind. molecules
    com_unwrap      = [side[0] for side in ctc_unwrap]
    traj_com_unwrap = [side[1] for side in ctc_unwrap]

    # current and previous atom positions in ind. leaflets
    mass_sys     = [side[0] for side in smpd_arrays]
    tmass_sys    = [side[1] for side in smpd_arrays]
    pos_sys      = [side[2] for side in smpd_arrays]
    pos_sys_prev = [side[3] for side in smpd_arrays]
    # displacements of atoms in individual leaflets
    displ_sys           = [side[4] for side in smpd_arrays]
    # COM displacement of individual leaflets
    displ_sys_com       = [side[5] for side in smpd_arrays]

    # unwrapped leaflet COM
    com_sys_unwrap      = [side[0] for side in ctc_sys_unwrap]
    # trajectory of unwrapped leaflet COM
    traj_com_sys_unwrap = [side[1] for side in ctc_sys_unwrap]

    # UNWRAPPING
    print('# UNWRAP trajectories')
    # sys.exit(0)
    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]  # Is this update frame ?

        # do frame-wise bilayer recentering - remaining translation
        if center:
            Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
            zref = np.zeros([2], dtype=float)
            for i in range(0, nside):
                zref[i] = np.mean(Lag_ref[i].positions[:, 2])
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
            pos_sys[j] = ag_sys[j].positions

            if i == 0:
                # Initialization

                # positions, unwrapped COM, and traj of unwrapped COM for the entire system
                mymsd.init_unwrap_sys_com(
                    pos_sys[j], mass_sys[j], tmass_sys[j], pos_sys_prev[j],
                    com_sys_unwrap[j], traj_com_sys_unwrap[j])

                # positions, unwrapped COM, and traj of unwraped COM for individual molecules
                mymsd.init_unwrap_mol_com(pos[j], mass_mol[j], tmass_mol[j], pos_prev[j],
                                          pos_unwrap[j], com_unwrap[j], traj_com_unwrap[j])
                # print(f'# leaflet {sside[j]}: init unwrapped com/traj done')
            else:
                # get sys COM displacement and update positions of system (in pos_sys_prev)
                # Current positions of system is obtained inside the function
                mymsd.calculate_displ_sys_com(i, box, pos_sys[j], pos_sys_prev[j],
                                              mass_sys[j], tmass_sys[j], displ_sys_com[j])
                # print(f'# leafelt {sside[j]}: sys COM displacement:',displ_sys_com[j])

                # update unwrapped system COM
                com_sys_unwrap[j] = com_sys_unwrap[j] + displ_sys_com[j]
                # update unwrapped system COM trajectory
                np.copyto(traj_com_sys_unwrap[j][i], com_sys_unwrap[j])

                # update unwrapped mol positions
                mymsd.update_unwrapped_mol_pos(i, box, pos[j], pos_prev[j],
                                               pos_unwrap[j], displ_sys_com[j])

                # calculate mol COMs & update their traj
                mymsd.update_unwrapped_mol_com_traj(i, pos_unwrap[j], mass_mol[j], tmass_mol[j],
                                                    com_unwrap[j], traj_com_unwrap[j])

                # print(com_sys_unwrap[j])
                # print(com_unwrap[j])

    print('# UNWRAPPING TRAJ & COMS DONE')
    # sys.exit(0)

    if qcomsys:
        print('# Write unwrapped COM of leaflets')
        write_com_sys(traj_com_sys_unwrap, nside,
                      framenum, interval, odir, sside)

    if qcommol:
        print('# Write unwrapped COMs of individual molecules')
        write_com_mol(traj_com_unwrap, nside, nmol,
                      framenum, interval, odir, sside)

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
            taus, framenum, traj_com_unwrap[i], id_type[i], msd[i])

    # Write MSD outputs
    if qb:
        # calculate MSD over bilayers
        bmsd = mymsd.calculate_msd_bilayer(msd, nside, ntype, ntau, nmol_type)
        write_msd_outputs_bilayer(bmsd, taus, ntype, name_type, odir)
    else:  # if (qb is False):
        write_msd_outputs_leaflet(
            msd, taus, nside, ntype, name_type, odir, sside)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_msd_membrane(**settings)


if __name__ == '__main__':
    main()
