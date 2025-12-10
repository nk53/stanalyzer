#!/usr/bin/python
import argparse
import re
import typing as t
from collections.abc import Sequence

import numpy as np
import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta
from . import mol_atomgroup_util as mymol
from . import msd_util as mymsd

if t.TYPE_CHECKING:
    import numpy.typing as npt

ANALYSIS_NAME = 'msd_solution'

NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]


def process_args(sel: str, split_to_mol: str | None) -> ProcessedArgs:
    """
    ----------
    Process arguements
    ----------
    """

    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()
    # print("selection",selection);print(ntype); sys.exit(0);

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
    else:  # get my_dict["split"] upto ntype
        qsplit = []
        for i in range(0, ntype):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)
    return ProcessedArgs(selection, ntype, qsplit)


# Write system COM
def write_sys_com(traj_com_sys_unwrap: NDFloat64,
                  framenum: int, interval: int, time_step: float,
                  odir: str, suffix: str) -> None:
    sout = '#     frame       COM\n'
    sout += '#                   x           y          z\n'
    for j in range(0, framenum):
        tcom = traj_com_sys_unwrap[j]
        sout += f' {time_step*(interval*j+1):10.5f}'
        # sout += f' {interval*j+1:10d}'
        sout += f' {tcom[0]:10.5f} {tcom[1]:10.5f} {tcom[2]:10.5f}\n'
    # print(sout)
    sta.write_to_outfile(f'{odir}/sys_com_{suffix}.dat', sout)


def write_com_mol(traj_com_unwrap: NDFloat64 | list[NDFloat64], nmol: int,
                  framenum: int, interval: int, time_step: float,
                  odir: str, suffix: str) -> None:
    """
    ----------
    Write unwrapped COMs of individual molecules
    ----------
    """

    sout = '#     frame       COM\n'
    sout += '#                   x           y          z\n'
    for i in range(0, framenum):
        tcom = traj_com_unwrap[i]
        sout += f' {time_step*(interval*i+1):10.5f}'
        # sout += f' {interval*i+1:10d}'
        for j in range(0, nmol):
            for k in range(0, 3):
                sout += f' {tcom[j, k]:10.5f}'
        sout += '\n'
    # print(sout)
    sta.write_to_outfile(f'{odir}/mol_com_{suffix}.dat', sout)


# Write x,y, and z-components of MSD for given molecule type
def write_msd(time_step: float,
              msd: NDFloat64, taus: list[int], fout: str) -> None:
    """
    ----------
    Write x,y, & z-components of MSD for a given molecule type"
    ----------
    """

    sout = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'
    ntau = len(taus)
    for i in range(0, ntau):
        sout += f'{time_step*taus[i]:10.5f}'
        # sout += f'{taus[i]:10d}'
        for j in range(0, 3):  # x,y,z
            sout += f' {msd[i, j]:10.5f}'
        sout += '\n'
    sta.write_to_outfile(fout, sout)


# Write MSD outputs
def write_msd_outputs(time_step: float,
                      msd: Sequence[NDFloat64] | NDFloat64, taus: list[int], ntype: int,
                      name_type: list[str], numb_type: list[int],
                      odir: str, suffix: str) -> None:
    for i in range(0, ntype):
        tname = name_type[i]
        print(f'# Write MSD output for {tname}')

        if numb_type[i] == 0:  # NA
            fout = f'{odir}/NA_{tname.lower()}_{suffix}.dat'
            sout = f'no {tname} in the system.\n'
            sta.write_to_outfile(fout, sout)
        else:
            fout = f'{odir}/{tname.lower()}_{suffix}.dat'
            write_msd(time_step,
                      msd[i], taus, fout)


def write_mol_info(nmol: int, id_type: list[int], name_type: list[str],
                   odir: str, suffix: str) -> None:
    # write molecule info: molecule number name
    sout = '# mol.index     name\n'
    # loop over molecules
    for j in range(0, nmol):
        jtype = id_type[j]  # molecule type index of molecule, j
        sout += f' {j:10d} {name_type[jtype]}\n'
    fout = f'{odir}/mol_info_{suffix}.dat'
    sta.write_to_outfile(fout, sout)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval', 'time_step')
    parser.add_argument(
        '--sel', metavar='selection',
        help='Selection for individual molecule. MSD will be calculated for COMs of individual '
             'molecules.')
    parser.add_argument(
        '--split',
        help='Y/N. If N, the atom group for the selection is considered as that for a single '
             'molecule. If Y, the atom group is further splitted to molecular level based on '
             'segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--qdrift', metavar='OPT', default='disabled',
                        choices=['enabled', 'disabled'],
                        help='If set on, COM drift is corrected.')
    parser.add_argument('--suffix', type=str, default='0',
                        help='Suffix  to output file(s)')
    parser.add_argument(
        '--qcomsys', default=False, action='store_true',
        help='If set to True, write unwrapped COM trajectory of system.')
    parser.add_argument(
        '--qcommol', default=False, action='store_true',
        help='If set to True, write unwrapped COM trajectoreis for molecules.')

    return parser


def run_msd_solution(sel: str, split: str | None, qdrift: str, qcomsys: bool, qcommol: bool,
                     psf: sta.FileRef, traj: sta.FileRefList, time_step: float | str,
                     suffix: str, interval: int = 1) -> None:
    """
    ----------
    Calculate Mean Square Displacement of COMs of selected molecule types
    ----------
    """

    # process non-system arguments
    selection, ntype, qsplit = process_args(sel, split)
    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    # print summary of arguments
    for i in range(0, ntype):
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    print('Suffix to output files', suffix)
    print('COM drift correction:', qdrift)
    print('Writing unwrapped system COM:', qcomsys)
    print('Writing unwrapped mol. COMs:', qcommol)
    print(f'MSD will be caulated every {interval} frames in delay time')

    # output dir
    odir = "./"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    # Generate molecule groups
    name_type, nmol_type, nmol, id_type, ag =\
        mymol.generate_mol_groups(u, ntype, selection, qsplit)

    # Set arrays in use
    # For molecules
    # Set mass, position, and displacement arrays
    # Will be input for calculation of unwrapped COM of individual molecules
    mass_mol, tmass_mol, pos, pos_prev, displ, pos_unwrap =\
        mymsd.set_mass_pos_displ_arrays(nmol, ag)

    # Set unwrapped COM positions and their trajectories for individual molecules
    com_unwrap, traj_com_unwrap =\
        mymsd.setup_unwrapped_mol_com_traj_array(ag, framenum)

    # For System
    if qdrift == 'enabled':
        # There should be no system COM drift
        ag_sys = u.atoms
        mass_sys, tmass_sys, pos_sys, pos_sys_prev, displ_sys, displ_sys_com = \
            mymsd.setup_sys_mass_pos_displ_arrays(ag_sys)

        com_sys_unwrap, traj_com_sys_unwrap = \
            mymsd.setup_unwrapped_com_traj_array(framenum)
    else:
        displ_sys_com = np.zeros([3], dtype=float)  # Set no COM drift

    # UNWRAPPING
    print('# UNWRAP trajectories')
    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]

        # get box size
        xtla, xtlb, xtlc = ts.dimensions[:3]
        box = np.array([xtla, xtlb, xtlc], dtype=float)

        # read cooridnates
        mymsd.read_coor(nmol, pos, ag)

        if i == 0:
            # Initialization
            if qdrift == 'enabled':
                # pos, unwrapped COM, and traj. of unwrapped COM for the system
                mymsd.init_unwrap_sys_com(pos_sys, mass_sys, tmass_sys, pos_sys_prev,
                                          com_sys_unwrap, traj_com_sys_unwrap)

            # positions, unwrapped COM, and traj of unwraped COM for individual molecules
            mymsd.init_unwrap_mol_com(pos, mass_mol, tmass_mol, pos_prev,
                                      pos_unwrap, com_unwrap, traj_com_unwrap)
            print('# init unwrap_pos done')
            # sys.exit(0)
        else:
            if qdrift == 'on':
                # calculate sys COM displacement & update pos_sys
                mymsd.calculate_displ_sys_com(i, box, pos_sys, pos_sys_prev,
                                              mass_sys, tmass_sys, displ_sys_com)
                # print(f'sys COM displacement:',displ_sys_com)

                # update unwrapped system COM
                com_sys_unwrap = com_sys_unwrap + displ_sys_com
                # update unwrapped system COM trajectory
                np.copyto(traj_com_sys_unwrap[i], com_sys_unwrap)

            # update unwrapped mol positions
            mymsd.update_unwrapped_mol_pos(
                i, box, pos, pos_prev, pos_unwrap, displ_sys_com)

            # calculate mol COMs & update their traj
            mymsd.update_unwrapped_mol_com_traj(
                i, pos_unwrap, mass_mol, tmass_mol, com_unwrap, traj_com_unwrap)

    print('# UNWRAPPING TRAJ & COMS DONE')
    # sys.exit(0)

    if qcomsys and qdrift == 'enabled':
        print('# Write unwrapped system COM')
        write_sys_com(traj_com_sys_unwrap, framenum, interval, time_step,
                      odir, suffix)

    if qcommol:
        print('# Write unwrapped COM of individual molecules')
        write_com_mol(traj_com_unwrap, nmol, framenum, interval, time_step,
                      odir, suffix)
        print('# Write mol. info : (mol. number,type)')
        write_mol_info(nmol, id_type, name_type, odir, suffix)

    print('# MSD calculations')
    taus = [int(interval*i) for i in range(0, framenum)]
    ntau = len(taus)  # number of data points along the delay time
    # DEBUG
    # print(f'ntau={ntau}')
    # print(taus)
    # for i in range(0,framenum):
    # print(f'i={interval*i}')
    # sys.exit(0)
    # ####

    # Setup msd arrays for individual molecule types
    msd = mymsd.setup_msd_arrays(ntype, ntau)

    # Calculate MSD for delay times, tau in {taus}
    mymsd.calculate_msd(taus, framenum, interval,
                        traj_com_unwrap, id_type, msd)

    # Write MSD outputs
    write_msd_outputs(time_step,
                      msd, taus, ntype, name_type, nmol_type, odir, suffix)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    # non-system arguments will be handled at the beginnig of this function
    run_msd_solution(**settings)


if __name__ == '__main__':
    main()
