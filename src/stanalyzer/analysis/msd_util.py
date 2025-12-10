#!/usr/bin/python
from collections.abc import Sequence
import typing as t

import numpy as np

if t.TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray
    from MDAnalysis import AtomGroup

T = t.TypeVar('T')
Tup2: t.TypeAlias = tuple[T, T]
NDFloat64: t.TypeAlias = 'NDArray[np.float64]'


class MassPosDisplTup1(t.NamedTuple):
    mass_mol:   list[NDFloat64]
    tmass_mol:  list[float]
    pos:        list[NDFloat64]
    pos_prev:   list[NDFloat64]
    displ:      list[NDFloat64]
    pos_unwrap: list[NDFloat64]


class MassPosDisplTup0(t.NamedTuple):
    mass_sys:       NDFloat64
    tmass_sys:      float
    pos_sys:        NDFloat64
    pos_sys_prev:   NDFloat64
    displ_sys:      NDFloat64
    displ_sys_com:  NDFloat64


def set_mass_pos_displ_arrays(nmol: int, ag: list['AtomGroup']) -> MassPosDisplTup1:
    """
    ----------
    Set mass, positions, and displacements of atoms in a atom group, ag.
    ----------
    input
          nmol      : number of molecules
          ag        : atom group arrays for individual molecules

    output
          mass_mol  : mass of ind. atoms in ind. molecules
          tmass_mol : mass of individual molecules
          pos       : current positions
          pos_prev  : previous positions
          displ     : displacements
          pos_unwrap: unwrapped positions
    """

    mass_mol: list[NDFloat64] = []
    tmass_mol: list[float] = []
    pos: list[NDFloat64] = []        # positions of each atom in each molecule
    # previous position of each atom in each molecule
    pos_prev: list[NDFloat64] = []
    # displacement array of each atoms in each molecule
    displ: list[NDFloat64] = []
    pos_unwrap: list[NDFloat64] = []  # unwrapped positions

    for i in range(0, nmol):
        natom = len(ag[i])
        # print(f'nmol = {nmol} imol={i} natom= {natom}');
        mass_mol.append(np.array([j.mass for j in ag[i]], dtype=float))
        tmass_mol.append(ag[i].total_mass(compound='group'))
        pos.append(np.zeros([natom, 3], dtype=float))
        pos_prev.append(np.zeros([natom, 3], dtype=float))
        displ.append(np.zeros([natom, 3], dtype=float))
        pos_unwrap.append(np.zeros([natom, 3], dtype=float))

    return MassPosDisplTup1(mass_mol, tmass_mol, pos, pos_prev, displ, pos_unwrap)


def setup_sys_mass_pos_displ_arrays(ag_sys: 'AtomGroup') -> MassPosDisplTup0:
    """
    ----------
    Setup systes mass,positions, and displacements and COM displacement
    ----------

    input
          ag_sys       : atom groups of system

    output
          mass_sys     : masses of atoms in the system
          tmass_sys    : total mass of the system
          pos_sys      : positions of atoms in the system
          pos_sys_prev : previous positions of atoms
          displ_sys    : displacements of atom positions
          displ_sys_com: displacement of COM position
    """
    natom = len(ag_sys)
    mass_sys = np.array([j.mass for j in ag_sys], dtype=float)
    tmass_sys = ag_sys.total_mass(compound='group')
    pos_sys = np.zeros([natom, 3], dtype=float)
    pos_sys_prev = np.zeros([natom, 3], dtype=float)
    displ_sys = np.zeros([natom, 3], dtype=float)
    displ_sys_com = np.zeros([3], dtype=float)

    return MassPosDisplTup0(mass_sys, tmass_sys, pos_sys, pos_sys_prev, displ_sys, displ_sys_com)


# Read coordinates of molecules
def read_coor(nmol: int, pos: list['NDArray'], ag: 'AtomGroup') -> None:
    """
    ----------
    Read positions of atom groups
    ----------

    input
          nmol: total number of molecules
          pos : positions
          ag  : atom groups

    output
          pos : updated positions
    """

    for i in range(0, nmol):
        tpos: 'ArrayLike' = ag[i].positions
        np.copyto(pos[i], tpos)


def setup_unwrapped_com_traj_array(framenum: int) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Setup arrays for unwrapped COM and trajectories
    ----------

    input
          framenum: number of frames in input trajectories

    output
          com     : COM
          traj_com: COM trajectory
    """
    com = np.zeros([3], dtype=float)
    traj_com = np.zeros([framenum, 3], dtype=float)
    return com, traj_com


def setup_unwrapped_mol_com_traj_array(
        ag: list['AtomGroup'], framenum: int) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Setup arrays for unwrapped positions and their trajectories
    ----------

    input
          ag      : atom groups
          framenum: number of frames in input trajectories

    output
          com_mol : COMs of individual molecules
          traj_com: trajectory of com_mol
    """

    nmol = len(ag)
    com_mol = np.zeros([nmol, 3], dtype=float)
    traj_com_mol = np.zeros([framenum, nmol, 3], dtype=float)
    return com_mol, traj_com_mol


def calculate_com(pos: 'NDArray', mass: 'ArrayLike', tmass: float) -> 'NDArray':
    """
    ----------
    Calculate COM of positions
    ----------

    input
          pos : positions
          mass: associated masses
          tmas: total mass

    output
          com : COM
    """

    tcom: 'NDArray' = (pos.T * mass).T   # numerator (vectors to be summed)
    com: 'NDArray' = np.sum(tcom, axis=0)/tmass

    return com


def init_unwrap_mol_com(pos: Sequence['NDArray'] | 'NDArray',
                        mass_mol: Sequence['ArrayLike'],
                        tmass_mol: NDFloat64 | Sequence[float],
                        pos_prev: Sequence['NDArray'],
                        pos_unwrap: Sequence['NDArray'],
                        com_unwrap: 'NDArray',
                        traj_com_unwrap: 'NDArray') -> None:
    """
    ----------
    Initialize unwrapped COMs and associated trajectories for individual molecules
    ----------

    input
          pos            : positions of atoms in individual molecules
          mass_mol       : masses of atoms in individual molecules
          tmass_mol      : masses of individual molecules

    input/output
          pos_prev       : previous positions of atoms in the system
          pos_unwrap     : unwrapped positions
          com_unwrap     : unwrapped system COM
          traj_com_unwrap: trajectory of unwrapped system COM
    """

    nmol = len(pos)
    for i in range(0, nmol):
        # initialize previous positions
        np.copyto(pos_prev[i], pos[i])
        np.copyto(pos_unwrap[i], pos[i])

        # calculate mol COM
        tcom = calculate_com(pos_unwrap[i], mass_mol[i], tmass_mol[i])
        # print(tcom)

        # initialize current unwrapped mol com & its trajectory at frame = 0
        np.copyto(com_unwrap[i, :], tcom)
        np.copyto(traj_com_unwrap[0, i, :], com_unwrap[i])

    # print(f'unwrapped frame 0: molid {nmol-1}',com_unwrap[-1])


def init_unwrap_sys_com(pos_sys: 'NDArray', mass_sys: 'ArrayLike',
                        tmass_sys: float, pos_sys_prev: 'NDArray',
                        com_sys_unwrap: 'NDArray',
                        traj_com_sys_unwrap: 'NDArray' | Sequence['NDArray']) -> None:
    """
    ----------
    Initialize unwrapped COM and associated trajectory for system
    ----------


    input
          pos_sys            : positions of atoms in the system
          mass_sys           : masses of atoms in the system
          tmass_sys          : total mass of the system

    input/output
          pos_sys_prev       : previous positions of atoms in the system
          com_sys_unwrap     : unwrapped system COM
          traj_com_sys_unwrap: trajectory of unwrapped system COM
    """

    # natom = len(pos_sys)  # for debugging?
    # initialize previous system pos
    np.copyto(pos_sys_prev, pos_sys)

    # calculate sys COM
    tcom = calculate_com(pos_sys, mass_sys, tmass_sys)

    # initialize current unwrapped sys com & associated trajectory at frame = 0
    np.copyto(com_sys_unwrap, tcom)
    np.copyto(traj_com_sys_unwrap[0], com_sys_unwrap)


def calculate_displ_sys_com(iframe: int, box: 'NDArray', pos_sys: 'NDArray',
                            pos_sys_prev: 'NDArray', mass_sys: 'NDArray',
                            tmass_sys: float, displ_sys_com: 'NDArray') -> None:
    """
    ----------
    Calculate COM displacement of the system
    ----------

    input
          iframe       : the current frame index
          box          : the current system sizes
          pos_sys      : positions of atoms in the system
          mass_sys     : masses of atoms in the system
          tmass_sys    : the total mass of the system

    input/output
          pos_sys_prev : previous positions of atoms in the system
          displ_sys_com: displacement of system COM
    """

    # natom = len(pos_sys)
    displ_sys = pos_sys - pos_sys_prev
    # tmpdispl = displ_sys; print(tmpdispl)
    displ_sys = displ_sys - np.sign(np.trunc(displ_sys/(box/2.0))) * box
    # tmpdispl = displ_sys; print(tmpdispl)

    # # DEBUG
    # maxdispl = np.max(np.absolute(displ_sys))  # max. displacement
    # if maxdispl > 10.0:
    #     print(f'frame {iframe}: max. displacement {maxdispl:10.5f} box:', box)
    #     # print(f'displacement:',displ)
    # # DEBUG

    # calculation of displacement of system COM
    displ_sys = (displ_sys.T * mass_sys).T  # numerator of disp_sys_com
    # print(displ_sys)
    tdispl_com = np.sum(displ_sys, axis=0)/tmass_sys
    np.copyto(displ_sys_com, tdispl_com)
    # print(displ_sys_com)

    # update previous positions
    np.copyto(pos_sys_prev, pos_sys)


def update_unwrapped_mol_pos(iframe: int, box: 'NDArray', pos: list['NDArray'],
                             pos_prev: list['NDArray'], pos_unwrap: list['NDArray'],
                             displ_sys_com: 'NDArray') -> None:
    """
    ----------
    Update unwrapped positions of individual molecules
    ----------
    NOTE: COM drift (displ_sys_com) is corrected


    input
          iframe       : the current frame index
          box          : the current system sizes
          pos          : current position
          displ_sys_com: displacement of system COM

    input/output
          pos_prev     : previous positions of atoms in individual molecules
          pos_unwrap   : unwrapped atom positions
    """

    nmol = len(pos)
    for i in range(0, nmol):
        displ = pos[i] - pos_prev[i]
        displ = displ - np.sign(np.trunc(displ/(box/2))) * box

        # # DEBUG
        # maxdispl = np.max(np.absolute(displ))
        # if maxdispl > 10.0:
        #     print(
        #         f'# frame {iframe}: imol= {i} max. displ. {maxdispl:10.5f} box:', box)
        #     print(displ)
        # # DEBUG

        # COM drift correction
        for j in range(0, len(displ)):
            displ[j] = displ[j] - displ_sys_com
        tpos_unwrap = pos_unwrap[i] + displ
        np.copyto(pos_unwrap[i], tpos_unwrap)

        # update previous position
        np.copyto(pos_prev[i], pos[i])

    # print(f'# displ:',displ)


def update_unwrapped_mol_com_traj(iframe: int,
                                  pos_unwrap: list['NDArray'],
                                  mass_mol: list['NDArray'],
                                  tmass_mol: list[float],
                                  com_unwrap: 'NDArray',
                                  traj_com_unwrap: 'NDArray') -> None:
    """
    ----------
    Update unwrapped COMs of individual molecules
    ----------

    input
          iframe          : the current frame index
          pos_unwrap      : current unwraped position
          mass_mol        : masses of atoms in individual molecules
          tmass_mol       : masses of individual molecules

    input/output
          com_unwrap      : unwrapped COMs of individual molecules
          traj_com_unwrap : trajectory of unwrapped COMs
    """

    nmol = len(pos_unwrap)
    # calculate COM of individual molecules
    for i in range(0, nmol):
        tcom = calculate_com(pos_unwrap[i], mass_mol[i], tmass_mol[i])
        np.copyto(com_unwrap[i], tcom)
    # print(f'unwrapped frame {iframe}: molid {nmol-1}',com_unwrap[-1])

    # update unwrapped trajectory & previous positions
    np.copyto(traj_com_unwrap[iframe, :, :], com_unwrap)


def setup_msd_arrays(ntype: int, ntau: int) -> NDFloat64:
    """
    ----------
    Setup MSD arrays for individual molecule types
    ----------

    input
          ntype: number of molecule types
          ntau : number of lag time points

    output
          msd  : time series of x,y,& z components of MSD for individual molecule types
    """

    msd = np.zeros([ntype, ntau, 3], dtype=float)  # msd[type,tau,x/y/z]

    return msd


def calculate_msd_tau(tau: int,
                      framenum: int, interval: int,
                      ntype: int,
                      id_type: Sequence[int],
                      traj_com_unwrap: NDFloat64) -> NDFloat64:
    """
    ----------
    Calculate MSD at a given lag time, tau
    ----------

    input
          tau            : lag time
          framenum       : ST-analyaer processed number of frames from input trajectories
          interval       : ST-analyzer set frame interval for analysis
          ntype          : number of molecule types
          id_type        : molecule type indices of individual molecules
          traj_com_unwrap: unwrapped COM trajectories of individual molecules

    output
          tmsd           : x,y, & z-components of MSD at tau for individual molecule types
    """

    # set msd array
    tmsd = np.zeros([ntype, 3], dtype=float)

    nmol = len(id_type)

    # set square displacement arrays
    sd_arr: list[list[list[NDFloat64]]] = []
    for i in range(0, ntype):
        sd_arr.append([])
        for j in range(0, 3):
            sd_arr[i].append([])

    # calculate square displacement of individual molecules
    # & update square displacement array
    ttau = int(tau/interval)
    for i in range(0, framenum - ttau):
        dis: NDFloat64 = traj_com_unwrap[i +
                                         ttau, :, :] - traj_com_unwrap[i, :, :]
        dis = np.square(dis)
        # loop over lipids & get msd for individual lipid type
        for j in range(0, nmol):
            k = id_type[j]  # molecule type index
            tsd: NDFloat64 = dis[j]
            for m in range(0, 3):
                sd_arr[k][m].append(tsd[m])

    # get MSD for individual molecule types
    for i in range(0, ntype):
        for j in range(0, 3):
            tmsd[i][j] = np.mean(sd_arr[i][j])

    return tmsd


def calculate_msd(taus: list[int], framenum: int, interval: int,
                  traj_com_unwrap: NDFloat64, id_type: Sequence[int],
                  msd: NDFloat64) -> None:
    """
    ----------
    Calculate MSD for given set of delay times
    ----------

    input
          taus: delay times
          framenum: ST-analyzer processed number of frames from input trajectories
          interval: ST-analyer set frameinterval for analysis
          traj_com_unwrap: unwrapped COM trajectories of individual molecules
          id_type        : molecule type indices of individual molecules

    input/output
          msd            : x,y, & z-components of MSD
    """

    ntau = len(taus)
    ntype = len(msd)
    for i in range(0, ntau):
        tau = taus[i]
        if int(tau/interval) > framenum - 1:
            break
        if tau % 100 == 0:
            print(f'MSD progress {tau}/{taus[-1]}')

        # Calculate MSD(tau)
        tmsd = calculate_msd_tau(tau, framenum, interval, ntype,
                                 id_type, traj_com_unwrap)

        # update MSD
        np.copyto(msd[:, i, :], tmsd)


def calculate_msd_bilayer(msd: Sequence[NDFloat64], nside: int, ntype: int,
                          ntaus: int, nmol_type: Sequence[Sequence[int]]) -> NDFloat64:
    """
    ----------
    Calculate bilayer msd from leaflet msd
    ----------

    input
          msd      : leaflet MSD, [nside][ntype,ntaus,3]
          nside    : number of leaflets, 2
          ntype    : number of unique molecule types
          ntaus    : number of lag time data points
          nmol_type: number of unique molecule types in each leafelt

    output
          bmsd     : MSD from both leaflets
    """
    bmsd = np.zeros([ntype, ntaus, 3], dtype=float)
    # sum over leaflets - nomralization factor
    bnmol_type = np.sum(nmol_type, axis=0)

    # -----------------------------------------
    # Calculate sum of square displacement
    # -----------------------------------------
    # leaflet data: msd[side][type,tau]
    #
    # sd[type,tau] =
    #   msd[up][type,tau] * nmol_type[up][type] + msd[dn][type,tau] * nmol_type[dn][type]
    #
    # msd[type,tau,:] = sd[type,tau]/np.sum(nmol_type,axis = 0)
    #
    for j in range(0, ntype):
        tnmolj = bnmol_type[j]
        if tnmolj == 0:
            continue  # when there's no type j in bilayer

        for k in range(0, ntaus):
            for i in range(0, nside):
                tnmolij = nmol_type[i][j]
                if tnmolij == 0:
                    continue  # when there's no type j in leaflet i

                for m in range(0, 3):
                    bmsd[j, k, m] += msd[i][j, k, m] * \
                        nmol_type[i][j] / bnmol_type[j]

    return bmsd
