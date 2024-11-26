import argparse
import typing as t

import MDAnalysis as mda    # type: ignore
import numpy as np
from MDAnalysis.analysis import distances  # type: ignore

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_int, p_float

ANALYSIS_NAME = 'contacts'

CPair: t.TypeAlias = tuple[str, int, str, int]
CDict: t.TypeAlias = dict[CPair, int]


def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Frequency"
    else:
        header_str = "Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Frequency"

    print(header_str, file=outfile)

    return header_str


def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   out: sta.FileRef, contact_threshold: float = 5.0,
                   debug: bool = False) -> None:
    """Writes contacts to `out` file"""

    if contact_threshold <= 0:
        raise ValueError("contact_threshold must be a positive number, "
                         f"not '{contact_threshold}'")

    universe = mda.Universe(psf, traj)

    # Create a dictionary to store contact frequency information
    contact_frequency: CDict = {}
    residues = universe.select_atoms(sel).residues

    # Iterate over frames in the trajectory
    for ts in universe.trajectory:
        # Select all atoms
        atoms = universe.select_atoms(sel)

        # Calculate the pairwise distances
        if debug:
            pairwise_distances = distances.distance_array(
                atoms.positions, atoms.positions, box=atoms.dimensions
            )
            print(pairwise_distances)

        # Check distances between residues
        num_residues = len(residues)
        for i in range(num_residues):
            for j in range(i + 1, num_residues):
                # Get the residue positions
                residue_i = residues[i].atoms
                residue_j = residues[j].atoms

                # Calculate distance between the center of mass of residues
                com_i = residue_i.center_of_mass()
                com_j = residue_j.center_of_mass()

                # Calculate distance
                dist = np.linalg.norm(com_i - com_j)

                if dist < contact_threshold:
                    contact = (residues[i].resname, residues[i].resid,
                               residues[j].resname, residues[j].resid)

                    # Increment the contact frequency
                    if contact in contact_frequency:
                        contact_frequency[contact] += 1
                    else:
                        contact_frequency[contact] = 1

    # Write the contact frequencies to the output file
    with sta.resolve_file(out, 'w') as outfile:
        outfile.write(
            "# Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Frequency\n")
        for contact, freq in contact_frequency.items():
            res_i_name, res_i_id, res_j_name, res_j_id = contact
            outfile.write(
                f"{res_i_name} {res_i_id} {res_j_name} {res_j_id} {freq}\n")

    print(f"Contact frequencies written to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('-rp', '--ref-psf', '--ref-psf-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('-rc', '--ref-coor', '--ref-coor-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="Coordinate file to use for reference, if not same as --traj")
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for contact calculation")
    parser.add_argument('--contact-threshold', type=p_float, metavar='N', default='5.0',
                        help="Distance cutoff for calculating the contact frequency.")
    parser.add_argument('-rt', '--ref-frame-type', choices=['specific', 'average'],
                        default='specific', metavar='TYPE',
                        help="specific: use a given frame as the reference; "
                        "average: use average structure")
    parser.add_argument('-rn', '--ref-frame-num', type=p_int, default=1, metavar='N',
                        help="Frame to use for reference coordinates (default: 1). "
                        "Only meaningful if --ref-frame-type is 'specific'")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()
