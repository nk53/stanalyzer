import argparse
import typing as t

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import self_capped_distance

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

ANALYSIS_NAME = 'contacts'

CPair: t.TypeAlias = tuple[str, int, str, int]
CDict: t.TypeAlias = dict[CPair, int]

def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    header_str = "Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Frequency"

    if not np_formatted:
        header_str = "# " + header_str

    if outfile is not None:
        print(header_str, file=outfile)

    return header_str

def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   out: sta.FileRef, contact_threshold: float = 5.0,
                   interval: int = 1, debug: bool = False) -> None:
    """
    Calculate residue-residue contact frequencies.

    A contact is recorded when the distance between the centers of mass
    of two residues is below `contact_threshold`.

    Optimized version:
    - Computes residue COMs once per frame.
    - Uses MDAnalysis self_capped_distance() to find only nearby pairs.
    - Avoids O(N^2) Python distance checks for every residue pair.
    """

    if contact_threshold <= 0:
        raise ValueError(
            "contact_threshold must be a positive number, "
            f"not '{contact_threshold}'"
        )

    universe = mda.Universe(psf, traj)

    # Perform atom selection only once.
    # Residue objects remain valid as trajectory frames advance.
    atoms = universe.select_atoms(sel)
    residues = atoms.residues

    if len(residues) == 0:
        raise ValueError(f"No residues found for selection: {sel}")

    contact_frequency: CDict = {}

    # Cache residue metadata once instead of repeatedly accessing
    # res.resname and res.resid inside the frame loop.
    residue_labels: list[tuple[str, int]] = [
        (res.resname, res.resid) for res in residues
    ]

    # Iterate through trajectory frames.
    for step_num, ts in enumerate(universe.trajectory):

        if step_num % interval:
            continue

        # Compute each residue center-of-mass exactly once per frame.
        # Old implementation recomputed COMs for every residue pair.
        coms = np.asarray(
            [res.atoms.center_of_mass() for res in residues],
            dtype=np.float64,
        )

        # Use MDAnalysis spatial search to find only residue pairs
        # within the contact cutoff.
        #
        # This replaces the nested:
        #
        #   for i:
        #       for j:
        #           distance(...)
        #
        # loop from the original implementation.
        pairs, _ = self_capped_distance(
            coms,
            max_cutoff=contact_threshold,
            box=ts.dimensions,
            return_distances=True,
        )

        # Prevent duplicate counting within a frame.
        frame_contacts: set[CPair] = set()

        for i, j in pairs:

            if i == j:
                continue

            # Ensure contact ordering is consistent.
            if i > j:
                i, j = j, i

            res_i_name, res_i_id = residue_labels[i]
            res_j_name, res_j_id = residue_labels[j]

            frame_contacts.add(
                (res_i_name, res_i_id,
                 res_j_name, res_j_id)
            )

        if debug:
            print(f"frame={step_num} contacts={len(frame_contacts)}")

        # Increment contact frequency once per frame.
        for contact in frame_contacts:
            contact_frequency[contact] = (
                contact_frequency.get(contact, 0) + 1
            )

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)

        for contact, freq in contact_frequency.items():
            print(*contact, freq, file=outfile)

    print(f"Contact frequencies written to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for contact calculation")
    parser.add_argument('--contact-threshold', type=p_float, metavar='N', default='5.0',
                        help="Distance cutoff for calculating the contact frequency.")
    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()