#!/usr/bin/python
from collections.abc import Sequence
import typing as t

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
from scipy import fft as scipy_fft

if t.TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray
    from MDAnalysis import AtomGroup

T = t.TypeVar('T')
Tup2: t.TypeAlias = tuple[T, T]
NDFloat64: t.TypeAlias = 'NDArray[np.float64]'
NDIntp: t.TypeAlias = 'NDArray[np.intp]'


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


class PackedMoleculeData(t.NamedTuple):
    """Molecule indices, masses, and flat arrays shared by the packed per-frame MSD pipeline."""

    atom_indices: NDIntp
    molecule_starts: NDIntp
    atom_masses: NDFloat64
    molecule_masses: NDFloat64
    coordinate_buffer: NDFloat64
    pos: NDFloat64
    pos_prev: NDFloat64
    pos_unwrap: NDFloat64
    displacement: NDFloat64
    minimum_image_scratch: NDFloat64
    weighted_pos: NDFloat64
    mass_sums: NDFloat64
    com_unwrap: NDFloat64
    traj_com_unwrap: NDFloat64


def minimum_image_displacement(
        displacement: NDFloat64,
        box: NDFloat64,
        out: 'NDFloat64 | None' = None,
        scratch: 'NDFloat64 | None' = None) -> NDFloat64:
    """Apply minimum-image PBC for orthorhombic or triclinic boxes."""
    box_array = np.asarray(box, dtype=float)
    if box_array.shape not in {(3,), (6,)}:
        raise ValueError(
            "box must contain either three lengths or six unit-cell values"
        )
    if not np.isfinite(box_array).all() or np.any(box_array[:3] <= 0):
        raise ValueError("box dimensions must be finite and positive")

    if (
        len(box_array) == 6
        and not np.allclose(box_array[3:], 90.0)
    ):
        minimized = minimize_vectors(displacement, box_array)
        if out is None:
            return t.cast(NDFloat64, minimized)
        np.copyto(out, minimized)
        return out

    lengths = box_array[:3]
    if out is None:
        out = np.array(displacement, dtype=float, copy=True)
    elif out is not displacement:
        np.copyto(out, displacement)
    if scratch is None:
        scratch = np.empty_like(out)
    np.divide(out, lengths / 2.0, out=scratch)
    np.trunc(scratch, out=scratch)
    np.sign(scratch, out=scratch)
    np.multiply(scratch, lengths, out=scratch)
    np.subtract(out, scratch, out=out)
    return out


def setup_packed_molecule_data(
        ag: list['AtomGroup'],
        framenum: int) -> PackedMoleculeData:
    """Precompute molecule indices and allocate contiguous working arrays."""
    molecule_sizes = np.asarray(
        [len(molecule) for molecule in ag],
        dtype=np.intp,
    )
    molecule_starts = np.empty(len(ag), dtype=np.intp)
    if len(ag):
        molecule_starts[0] = 0
        if len(ag) > 1:
            np.cumsum(
                molecule_sizes[:-1],
                out=molecule_starts[1:],
            )
        atom_indices = np.concatenate(
            [molecule.indices for molecule in ag]
        ).astype(np.intp, copy=False)
        atom_masses = np.concatenate(
            [np.asarray(molecule.masses, dtype=float) for molecule in ag]
        )
        molecule_masses = np.asarray(
            [
                molecule.total_mass(compound='group')
                for molecule in ag
            ],
            dtype=float,
        )
        if (
            not np.isfinite(atom_masses).all()
            or np.any(atom_masses <= 0)
            or not np.isfinite(molecule_masses).all()
            or np.any(molecule_masses <= 0)
        ):
            raise ValueError(
                "All analyzed atoms and molecules require finite, "
                "positive masses"
            )
    else:
        atom_indices = np.empty(0, dtype=np.intp)
        atom_masses = np.empty(0, dtype=float)
        molecule_masses = np.empty(0, dtype=float)

    natoms = len(atom_indices)
    nmol = len(ag)
    coordinate_dtype = ag[0].positions.dtype if len(ag) else float
    return PackedMoleculeData(
        atom_indices=atom_indices,
        molecule_starts=molecule_starts,
        atom_masses=atom_masses,
        molecule_masses=molecule_masses,
        coordinate_buffer=np.empty((natoms, 3), dtype=coordinate_dtype),
        pos=np.empty((natoms, 3), dtype=float),
        pos_prev=np.empty((natoms, 3), dtype=float),
        pos_unwrap=np.empty((natoms, 3), dtype=float),
        displacement=np.empty((natoms, 3), dtype=float),
        minimum_image_scratch=np.empty((natoms, 3), dtype=float),
        weighted_pos=np.empty((natoms, 3), dtype=float),
        mass_sums=np.empty((nmol, 3), dtype=float),
        com_unwrap=np.empty((nmol, 3), dtype=float),
        traj_com_unwrap=np.empty((framenum, nmol, 3), dtype=float),
    )


def read_packed_coordinates(
        frame_positions: NDFloat64,
        data: PackedMoleculeData) -> None:
    """Gather all selected molecule coordinates with one indexed operation."""
    np.take(
        frame_positions,
        data.atom_indices,
        axis=0,
        out=data.coordinate_buffer,
    )
    np.copyto(data.pos, data.coordinate_buffer)


def calculate_packed_com(data: PackedMoleculeData) -> NDFloat64:
    """Calculate packed molecule COMs in legacy summation order.

    Returns a view of data.com_unwrap.
    """
    if not len(data.molecule_starts):
        return np.empty((0, 3), dtype=float)

    np.multiply(
        data.pos_unwrap,
        data.atom_masses[:, np.newaxis],
        out=data.weighted_pos,
    )
    np.add.reduceat(
        data.weighted_pos,
        data.molecule_starts,
        axis=0,
        out=data.mass_sums,
    )
    np.divide(
        data.mass_sums,
        data.molecule_masses[:, np.newaxis],
        out=data.com_unwrap,
    )
    return data.com_unwrap


def init_packed_molecule_com(data: PackedMoleculeData) -> None:
    """Initialize packed unwrapped coordinates and first-frame COMs."""
    np.copyto(data.pos_prev, data.pos)
    np.copyto(data.pos_unwrap, data.pos)
    np.copyto(data.com_unwrap, calculate_packed_com(data))
    np.copyto(data.traj_com_unwrap[0], data.com_unwrap)


def update_packed_molecule_com(
        iframe: int,
        box: NDFloat64,
        data: PackedMoleculeData,
        displ_sys_com: NDFloat64) -> None:
    """Update all unwrapped molecule positions and COMs for one frame."""
    np.subtract(data.pos, data.pos_prev, out=data.displacement)
    minimum_image_displacement(
        data.displacement,
        box,
        out=data.displacement,
        scratch=data.minimum_image_scratch,
    )
    np.subtract(
        data.displacement,
        displ_sys_com,
        out=data.displacement,
    )

    np.add(
        data.pos_unwrap,
        data.displacement,
        out=data.pos_unwrap,
    )
    np.copyto(data.pos_prev, data.pos)
    np.copyto(data.com_unwrap, calculate_packed_com(data))
    np.copyto(
        data.traj_com_unwrap[iframe],
        data.com_unwrap,
    )


def set_mass_pos_displ_arrays(nmol: int, ag: list['AtomGroup']) -> MassPosDisplTup1:
    """Allocate per-molecule mass, position, displacement, and unwrap arrays.

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
    pos: list[NDFloat64] = []
    pos_prev: list[NDFloat64] = []
    displ: list[NDFloat64] = []
    pos_unwrap: list[NDFloat64] = []

    for i in range(0, nmol):
        natom = len(ag[i])
        mass_mol.append(np.array([j.mass for j in ag[i]], dtype=float))
        tmass_mol.append(ag[i].total_mass(compound='group'))
        pos.append(np.zeros([natom, 3], dtype=float))
        pos_prev.append(np.zeros([natom, 3], dtype=float))
        displ.append(np.zeros([natom, 3], dtype=float))
        pos_unwrap.append(np.zeros([natom, 3], dtype=float))

    return MassPosDisplTup1(mass_mol, tmass_mol, pos, pos_prev, displ, pos_unwrap)


def setup_sys_mass_pos_displ_arrays(ag_sys: 'AtomGroup') -> MassPosDisplTup0:
    """Allocate system-wide mass, position, displacement, and COM-displacement 
arrays.

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


def read_coor(nmol: int, pos: list['NDArray'], ag: 'AtomGroup') -> None:
    """Copy per-molecule atom positions into preallocated arrays.

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
    """Allocate system COM and per-frame COM trajectory arrays.

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
    """Allocate per-molecule COM and per-frame COM-trajectory arrays.

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
    """Weighted center-of-mass from positions, masses, and total mass.

    input
          pos : positions
          mass: associated masses
          tmas: total mass

    output
          com : COM
    """

    tcom: 'NDArray' = (pos.T * mass).T
    com: 'NDArray' = np.sum(tcom, axis=0)/tmass

    return com


def init_unwrap_mol_com(pos: 'Sequence[NDArray] | NDArray',
                        mass_mol: Sequence['ArrayLike'],
                        tmass_mol: 'NDFloat64 | Sequence[float]',
                        pos_prev: Sequence['NDArray'],
                        pos_unwrap: Sequence['NDArray'],
                        com_unwrap: 'NDArray',
                        traj_com_unwrap: 'NDArray') -> None:
    """Seed per-molecule unwrapped positions, COMs, and trajectory at frame 0.

    input
          pos            : positions of atoms in individual molecules
          mass_mol       : masses of atoms in individual molecules
          tmass_mol      : masses of individual molecules

    input/output
          pos_prev       : previous positions of atoms in the system
          pos_unwrap     : unwrapped positions
          com_unwrap     : unwrapped COM of individual molecules
          traj_com_unwrap: trajectory of unwrapped COMs of individual molecules
    """

    nmol = len(pos)
    for i in range(0, nmol):
        np.copyto(pos_prev[i], pos[i])
        np.copyto(pos_unwrap[i], pos[i])

        tcom = calculate_com(pos_unwrap[i], mass_mol[i], tmass_mol[i])

        np.copyto(com_unwrap[i, :], tcom)
        np.copyto(traj_com_unwrap[0, i, :], com_unwrap[i])


def init_unwrap_sys_com(pos_sys: 'NDArray', mass_sys: 'ArrayLike',
                        tmass_sys: float, pos_sys_prev: 'NDArray',
                        com_sys_unwrap: 'NDArray',
                        traj_com_sys_unwrap: 'NDArray | Sequence[NDArray]') -> None:
    """Seed system unwrapped COM and trajectory at frame 0.

    input
          pos_sys            : positions of atoms in the system
          mass_sys           : masses of atoms in the system
          tmass_sys          : total mass of the system

    input/output
          pos_sys_prev       : previous positions of atoms in the system
          com_sys_unwrap     : unwrapped system COM
          traj_com_sys_unwrap: trajectory of unwrapped system COM
    """

    np.copyto(pos_sys_prev, pos_sys)
    tcom = calculate_com(pos_sys, mass_sys, tmass_sys)
    np.copyto(com_sys_unwrap, tcom)
    np.copyto(traj_com_sys_unwrap[0], com_sys_unwrap)


def calculate_displ_sys_com(iframe: int, box: 'NDArray', pos_sys: 'NDArray',
                            pos_sys_prev: 'NDArray', mass_sys: 'NDArray',
                            tmass_sys: float, displ_sys_com: 'NDArray',
                            scratch: 'NDArray | None' = None) -> None:
    """Compute mass-weighted system COM displacement for one frame.

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

    if scratch is None:
        displ_sys = pos_sys - pos_sys_prev
    else:
        displ_sys = scratch
        np.subtract(pos_sys, pos_sys_prev, out=displ_sys)
    displ_sys = minimum_image_displacement(
        displ_sys,
        box,
        out=displ_sys,
    )

    np.multiply(displ_sys, mass_sys[:, np.newaxis], out=displ_sys)
    tdispl_com = np.sum(displ_sys, axis=0)/tmass_sys
    np.copyto(displ_sys_com, tdispl_com)
    np.copyto(pos_sys_prev, pos_sys)


def update_unwrapped_mol_pos(iframe: int, box: 'NDArray', pos: list['NDArray'],
                             pos_prev: list['NDArray'], pos_unwrap: list['NDArray'],
                             displ_sys_com: 'NDArray') -> None:
    """Advance per-molecule unwrapped positions, correcting for system COM drift.

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
        displ = minimum_image_displacement(
            pos[i] - pos_prev[i],
            box,
        )

        # Apply the three-component COM drift to every atom at once.
        displ -= displ_sys_com
        tpos_unwrap = pos_unwrap[i] + displ
        np.copyto(pos_unwrap[i], tpos_unwrap)

        np.copyto(pos_prev[i], pos[i])


def update_unwrapped_mol_com_traj(iframe: int,
                                  pos_unwrap: list['NDArray'],
                                  mass_mol: list['NDArray'],
                                  tmass_mol: list[float],
                                  com_unwrap: 'NDArray',
                                  traj_com_unwrap: 'NDArray') -> None:
    """Recompute per-molecule COMs from unwrapped positions and record trajectory.

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
    for i in range(0, nmol):
        tcom = calculate_com(pos_unwrap[i], mass_mol[i], tmass_mol[i])
        np.copyto(com_unwrap[i], tcom)

    np.copyto(traj_com_unwrap[iframe, :, :], com_unwrap)


def setup_msd_arrays(ntype: int, ntau: int) -> NDFloat64:
    """Allocate per-type, per-lag x/y/z MSD array.

    input
          ntype: number of molecule types
          ntau : number of lag time points

    output
          msd  : time series of x,y,& z components of MSD for individual molecule types
    """

    msd = np.zeros([ntype, ntau, 3], dtype=float)

    return msd


def calculate_msd_tau(tau: int,
                      framenum: int, interval: int,
                      ntype: int,
                      id_type: Sequence[int],
                      traj_com_unwrap: NDFloat64) -> NDFloat64:
    """Compute per-type MSD at one lag time from unwrapped COM trajectories.

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

    type_ids = np.asarray(id_type, dtype=np.intp)
    type_indices = [
        np.flatnonzero(type_ids == molecule_type)
        for molecule_type in range(ntype)
    ]
    return _calculate_msd_tau(
        tau=tau,
        framenum=framenum,
        interval=interval,
        type_indices=type_indices,
        traj_com_unwrap=traj_com_unwrap,
    )


def _calculate_msd_tau(
        tau: int,
        framenum: int,
        interval: int,
        type_indices: Sequence[NDIntp],
        traj_com_unwrap: NDFloat64,
        scratch: 'NDFloat64 | None' = None) -> NDFloat64:
    """Calculate one lag-time MSD using vectorized time origins."""
    frame_lag = int(tau / interval)
    stop = framenum - frame_lag

    if scratch is None:
        squared_displacement = (
            traj_com_unwrap[frame_lag:framenum]
            - traj_com_unwrap[:stop]
        )
    else:
        squared_displacement = scratch[:stop]
        np.subtract(
            traj_com_unwrap[frame_lag:framenum],
            traj_com_unwrap[:stop],
            out=squared_displacement,
        )
    np.square(squared_displacement, out=squared_displacement)

    tmsd = np.empty((len(type_indices), 3), dtype=float)
    for molecule_type, indices in enumerate(type_indices):
        tmsd[molecule_type] = np.mean(
            squared_displacement[:, indices, :],
            axis=(0, 1),
        )

    return tmsd


def calculate_msd(taus: list[int], framenum: int, interval: int,
                  traj_com_unwrap: NDFloat64, id_type: Sequence[int],
                  msd: NDFloat64) -> None:
    """Fill per-type MSD for every requested lag using the direct engine.

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
    type_ids = np.asarray(id_type, dtype=np.intp)
    if traj_com_unwrap.shape[1] != len(type_ids):
        raise ValueError(
            "id_type length must match the molecule dimension of "
            "traj_com_unwrap"
        )
    if np.any(type_ids < 0) or np.any(type_ids >= ntype):
        raise ValueError("id_type contains an invalid molecule type index")

    type_indices = [
        np.flatnonzero(type_ids == molecule_type)
        for molecule_type in range(ntype)
    ]
    # Reuse one full-sized work array across all lags instead of allocating per lag.
    scratch = np.empty_like(traj_com_unwrap)

    for i in range(0, ntau):
        tau = taus[i]
        if int(tau/interval) > framenum - 1:
            break
        if tau % 100 == 0:
            print(f'MSD progress {tau}/{taus[-1]}')
        if tau == 0:
            msd[:, i, :] = 0.0
            continue

        tmsd = _calculate_msd_tau(
            tau=tau,
            framenum=framenum,
            interval=interval,
            type_indices=type_indices,
            traj_com_unwrap=traj_com_unwrap,
            scratch=scratch,
        )

        np.copyto(msd[:, i, :], tmsd)


def _next_fast_fft_length(minimum_length: int) -> int:
    """Return an efficient real-FFT length at least ``minimum_length``."""
    if minimum_length < 1:
        raise ValueError("minimum_length must be positive")
    return int(scipy_fft.next_fast_len(minimum_length, real=True))


def select_msd_engine(
        framenum: int,
        frame_lags: Sequence[int]) -> t.Literal['direct', 'fft']:
    """Select an engine from trajectory and requested-lag work estimates."""
    if framenum < 1:
        raise ValueError("framenum must be positive")
    lags = np.asarray(frame_lags, dtype=np.intp)
    if np.any(lags < 0) or np.any(lags >= framenum):
        raise ValueError("frame_lags contain a lag outside the trajectory")

    # Direct work scales with displacement samples; FFT work covers the full padded trajectory.
    direct_samples = int(np.sum(framenum - lags, dtype=np.int64))
    fft_length = _next_fast_fft_length(2 * framenum - 1)
    fft_work = fft_length * np.log2(fft_length)
    return 'fft' if direct_samples > fft_work else 'direct'


def calculate_msd_fft(
        taus: list[int], framenum: int, interval: int,
        traj_com_unwrap: NDFloat64, id_type: Sequence[int],
        msd: NDFloat64,
        chunk_size: int | None = None,
        workers: int = 1) -> None:
    """FFT autocorrelation MSD; faster for long trajectories but float-reduction order differs from the direct engine."""
    ntype = len(msd)
    type_ids = np.asarray(id_type, dtype=np.intp)
    if traj_com_unwrap.shape != (framenum, len(type_ids), 3):
        raise ValueError(
            "traj_com_unwrap shape must be (framenum, len(id_type), 3)"
        )
    if np.any(type_ids < 0) or np.any(type_ids >= ntype):
        raise ValueError("id_type contains an invalid molecule type index")

    frame_lags = np.asarray(taus, dtype=np.intp) // interval
    if np.any(frame_lags < 0) or np.any(frame_lags >= framenum):
        raise ValueError("taus contain a lag outside the trajectory")
    if chunk_size is not None and chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    if workers < 1:
        raise ValueError("workers must be positive")

    # Linear autocorrelation of N samples needs 2N-1 points; use a compact mixed-radix FFT length.
    minimum_fft_length = 2 * framenum - 1
    fft_length = _next_fast_fft_length(minimum_fft_length)
    for molecule_type in range(ntype):
        indices = np.flatnonzero(type_ids == molecule_type)
        if not len(indices):
            msd[molecule_type] = np.nan
            continue

        effective_chunk_size = (
            len(indices)
            if chunk_size is None
            else chunk_size
        )
        autocorrelation = np.zeros((framenum, 3), dtype=float)
        squared = np.zeros((framenum, 3), dtype=float)
        for chunk_start in range(0, len(indices), effective_chunk_size):
            chunk_indices = indices[
                chunk_start:chunk_start + effective_chunk_size
            ]
            coordinates = traj_com_unwrap[:, chunk_indices, :]
            spectrum = scipy_fft.rfft(
                coordinates,
                n=fft_length,
                axis=0,
                workers=workers,
            )

            # Convert F to |F|^2 in place to avoid a second complex FFT-sized array.
            np.square(spectrum.real, out=spectrum.real)
            np.square(spectrum.imag, out=spectrum.imag)
            spectrum.real += spectrum.imag
            spectrum.imag.fill(0.0)
            autocorrelation += scipy_fft.irfft(
                spectrum,
                n=fft_length,
                axis=0,
                workers=workers,
            )[:framenum].sum(axis=1)
            squared += np.einsum(
                'tmc,tmc->tc',
                coordinates,
                coordinates,
                optimize=True,
            )

        prefix = np.vstack(
            [np.zeros((1, 3), dtype=float), np.cumsum(squared, axis=0)]
        )
        for output_index, frame_lag in enumerate(frame_lags):
            origin_count = framenum - frame_lag
            numerator = (
                prefix[origin_count]
                + prefix[framenum]
                - prefix[frame_lag]
                - 2.0 * autocorrelation[frame_lag]
            )
            values = numerator / (origin_count * len(indices))
            # Roundoff in the FFT can make exact zeros very slightly negative.
            msd[molecule_type, output_index] = np.maximum(values, 0.0)

    msd[:, 0, :] = 0.0

# FFT MSD: all lags via C = IFFT(|FFT(r)|^2) (Fourier correlation theorem), O(N log N) vs O(N^2); zero-pad for linear autocorrelation; float roundoff differs from direct.

def calculate_msd_bilayer(msd: Sequence[NDFloat64], nside: int, ntype: int,
                          ntaus: int, nmol_type: Sequence[Sequence[int]]) -> NDFloat64:
    """Combine leaflet MSDs into a molecule-count-weighted bilayer MSD.

    input
          msd      : leaflet MSD, [nside][ntype,ntaus,3]
          nside    : number of leaflets, 2
          ntype    : number of unique molecule types
          ntaus    : number of lag time data points
          nmol_type: number of unique molecule types in each leafelt

    output
          bmsd     : MSD from both leaflets
    """
    msd_array = np.asarray(msd, dtype=float)
    counts = np.asarray(nmol_type, dtype=float)
    if msd_array.shape != (nside, ntype, ntaus, 3):
        raise ValueError("msd has an incompatible shape")
    if counts.shape != (nside, ntype):
        raise ValueError("nmol_type has an incompatible shape")

    totals = counts.sum(axis=0)
    weights = np.divide(
        counts,
        totals[np.newaxis, :],
        out=np.zeros_like(counts),
        where=totals[np.newaxis, :] != 0,
    )
    result = np.einsum('stkc,st->tkc', msd_array, weights, optimize=True)
    return t.cast(NDFloat64, result)
