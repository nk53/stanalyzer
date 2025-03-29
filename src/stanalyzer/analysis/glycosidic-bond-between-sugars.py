import argparse
from typing import Optional

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis import distances  # type: ignore[import-untyped]

ANALYSIS_NAME = 'glycosidic_bond_between_sugars'


def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time contacts"
    else:
        header_str = "#time contacts"

    print(header_str, file=outfile)

    return header_str


def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   out: sta.FileRef, ref_psf: Optional[sta.FileRef] = None,
                   ref_coor: sta.FileRef | None = None,
                   ref_frame_type: str = 'specific',
                   ref_frame_num: int = 1, interval: int = 1) -> None:
    """Writes contacts to `out` file"""

    # Load your trajectory and topology
    u = mda.Universe(psf, traj)

    # Define the selection strings for sugars and specific atoms
    # These should be adjusted based on specific system and atom names
    anomeric_carbon = u.select_atoms("  name C1")  # Example selection
    hydroxyl_group = u.select_atoms("  name O1")   # Example selection

    # Create an empty list to store bond distances
    distances_list = []

    # Iterate through frames in the trajectory
    for ts in u.trajectory:
        # Calculate the distance between the anomeric carbon and hydroxyl group
        dist = distances.distance_array(anomeric_carbon.positions, hydroxyl_group.positions)

        # Store the distance
        distances_list.append(dist)

    # Convert distances_list to a numpy array for easier analysis
    distances_array = np.array(distances_list)
    print(distances_array)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'pdb', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for glycosidic_bond calculation")
    parser.add_argument('-c', '--center', action='store_true')
    parser.add_argument('--sel', metavar='selection')

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()
