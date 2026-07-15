import typing as t

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import self_capped_distance

CPair: t.TypeAlias = tuple[str, int, str, int]
CDict: t.TypeAlias = dict[CPair, int]


def contacts_worker(task):
    (
        psf,
        traj,
        sel,
        contact_threshold,
        interval,
        debug,
        start,
        stop,
    ) = task

    universe = mda.Universe(psf, traj)
    atoms = universe.select_atoms(sel)
    residues = atoms.residues

    if len(residues) == 0:
        raise ValueError(f"No residues found for selection: {sel}")

    residue_atom_indices = [
        res.atoms.indices
        for res in residues
    ]

    all_masses = universe.atoms.masses

    residue_labels: list[tuple[str, int]] = [
        (res.resname, res.resid)
        for res in residues
    ]

    contact_frequency: CDict = {}

    for frame_idx in range(start, stop):
        if frame_idx % interval:
            continue

        ts = universe.trajectory[frame_idx]
        positions = universe.atoms.positions

        coms = np.empty(
            (len(residue_atom_indices), 3),
            dtype=np.float64,
        )

        for i, atom_idx in enumerate(residue_atom_indices):
            masses = all_masses[atom_idx]
            coords = positions[atom_idx]
            total_mass = masses.sum()

            coms[i] = (
                coords * masses[:, None]
            ).sum(axis=0) / total_mass

        pairs, _ = self_capped_distance(
            coms,
            max_cutoff=contact_threshold,
            box=ts.dimensions,
            return_distances=True,
        )

        frame_contacts: set[CPair] = set()

        for i, j in pairs:
            if i == j:
                continue

            if i > j:
                i, j = j, i

            res_i_name, res_i_id = residue_labels[i]
            res_j_name, res_j_id = residue_labels[j]

            frame_contacts.add(
                (
                    res_i_name,
                    res_i_id,
                    res_j_name,
                    res_j_id,
                )
            )

        if debug:
            print(f"frame={frame_idx} contacts={len(frame_contacts)}")

        for contact in frame_contacts:
            contact_frequency[contact] = (
                contact_frequency.get(contact, 0) + 1
            )

    return contact_frequency
