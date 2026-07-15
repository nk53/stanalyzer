from collections import defaultdict

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import normal, angle


PiPair = tuple[str, str]
PiEvents = dict[PiPair, set[int]]


RESIDUE_TO_RING_ATOMS = {
    "PHE": "name CG CD* CE* CZ",
    "TYR": "name CG CD* CE* CZ",
    "TRP": "name CD2 CE2 CZ2 CH2 CZ3 CE3",
    "HIS": "name CG ND1 CE1 NE2 CD2",
    "HSD": "name CG ND1 CE1 NE2 CD2",
    "HSE": "name CG ND1 CE1 NE2 CD2",
    "HSP": "name CG ND1 CE1 NE2 CD2",
    "HID": "name CG ND1 CE1 NE2 CD2",
    "HIE": "name CG ND1 CE1 NE2 CD2",
    "HIP": "name CG ND1 CE1 NE2 CD2",
    "CHID": "name CG ND1 CE1 NE2 CD2",
    "CHIE": "name CG ND1 CE1 NE2 CD2",
    "CHIP": "name CG ND1 CE1 NE2 CD2",
    "NHID": "name CG ND1 CE1 NE2 CD2",
    "NHIE": "name CG ND1 CE1 NE2 CD2",
    "NHIP": "name CG ND1 CE1 NE2 CD2",
}


CATION_SELECTION = (
    "(resname ARG LYS CARG CLYS NARG NLYS) "
    "and (name NE NH* NZ)"
)


def residue_label(atom) -> str:
    return f"{atom.segid}_{atom.resname}_{atom.resid}"


def pi_stacking_worker(task) -> PiEvents:
    (
        psf,
        traj,
        sel,
        pi_pi_dist_cutoff,
        pi_cation_dist_cutoff,
        interval,
        start,
        stop,
        debug,
    ) = task

    universe = mda.Universe(psf, traj)
    all_atoms = universe.select_atoms(sel)

    if len(all_atoms) == 0:
        raise ValueError(f"No atoms found for selection: {sel}")

    pi_rings = []
    pi_labels = []

    for residue in all_atoms.residues:
        atom_selection = RESIDUE_TO_RING_ATOMS.get(residue.resname)

        if atom_selection is None:
            continue

        ring = residue.atoms.select_atoms(atom_selection)

        if len(ring) < 3:
            continue

        pi_rings.append(ring)
        pi_labels.append(residue_label(ring[0]))

    if not pi_rings:
        raise ValueError("Unable to find aromatic residues")

    cations = all_atoms.select_atoms(CATION_SELECTION)

    cation_labels = [
        residue_label(atom)
        for atom in cations
    ]

    pi_cation_limit = np.pi / 6.0
    pi_cation_limit_opposite = np.pi - pi_cation_limit

    events: dict[PiPair, set[int]] = defaultdict(set)

    first_frame = start + (-start % interval)

    for frame_idx in range(first_frame, stop, interval):
        ts = universe.trajectory[frame_idx]

        ring_centers = np.empty(
            (len(pi_rings), 3),
            dtype=np.float64,
        )

        ring_normals = np.empty(
            (len(pi_rings), 3),
            dtype=np.float64,
        )

        for idx, ring in enumerate(pi_rings):
            positions = ring.positions
            center = positions.mean(axis=0)

            ring_centers[idx] = center

            v1 = positions[0] - center
            v2 = positions[1] - center
            ring_normals[idx] = normal(v1, v2)

        frame_events: set[PiPair] = set()

        # π–π candidates
        pi_pairs = capped_distance(
            ring_centers,
            ring_centers,
            max_cutoff=pi_pi_dist_cutoff,
            min_cutoff=1.0,
            box=ts.dimensions,
            return_distances=False,
        )

        for idx1, idx2 in pi_pairs:
            idx1 = int(idx1)
            idx2 = int(idx2)

            if idx1 >= idx2:
                continue

            frame_events.add(
                tuple(sorted(
                    (
                        pi_labels[idx1],
                        pi_labels[idx2],
                    )
                ))
            )

        # π–cation candidates
        if len(cations) > 0:
            pi_cation_pairs = capped_distance(
                ring_centers,
                cations.positions,
                max_cutoff=pi_cation_dist_cutoff,
                min_cutoff=1.0,
                box=ts.dimensions,
                return_distances=False,
            )

            cation_positions = cations.positions

            for ring_idx, cation_idx in pi_cation_pairs:
                ring_idx = int(ring_idx)
                cation_idx = int(cation_idx)

                center_to_cation = (
                    cation_positions[cation_idx]
                    - ring_centers[ring_idx]
                )

                radian = angle(
                    ring_normals[ring_idx],
                    center_to_cation,
                )

                aligned = (
                    -pi_cation_limit <= radian <= pi_cation_limit
                    or radian >= pi_cation_limit_opposite
                    or radian <= -pi_cation_limit_opposite
                )

                if aligned:
                    frame_events.add(
                        (
                            pi_labels[ring_idx],
                            cation_labels[cation_idx],
                        )
                    )

        for event in frame_events:
            events[event].add(frame_idx)

        if debug:
            print(
                f"frame={frame_idx} "
                f"pi_pi_candidates={len(pi_pairs)} "
                f"events={len(frame_events)}"
            )

    return dict(events)