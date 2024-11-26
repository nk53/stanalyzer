import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

ANALYSIS_NAME = 'helix_distance_crossing_angle'

def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time Distance Crossing_angle"
    else:
        header_str = "#time Distance Crossing_angle"

    print(header_str, file=outfile)

    return header_str


def write_helix_distance_crossing_angle(psf: sta.FileRef, traj: sta.FileRefList,
                                  helix1_start: int, helix1_end: int, 
                                  helix2_start: int, helix2_end: int, 
                                  out: sta.FileRef, interval: int = 1) -> None:
    """Writes distance and crossing angle between two helices to `out` file"""

    # Load the universe with the trajectory and topology
    mobile = mda.Universe(psf, traj)
    
    # Select atoms for helix1 and helix2 based on residue ranges
    helix1 = mobile.select_atoms(f"resid {helix1_start}:{helix1_end} and name CA")
    helix2 = mobile.select_atoms(f"resid {helix2_start}:{helix2_end} and name CA")

    # Function to calculate the principal axis
    def calculate_principal_axis(atoms):
        positions = atoms.positions - atoms.center_of_geometry()
        inertia_tensor = np.dot(positions.T, positions)
        eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
        principal_axis = eigenvectors[:, np.argmin(eigenvalues)]
        return principal_axis / np.linalg.norm(principal_axis)

    results = []
    for ts in mobile.trajectory[::interval]:
        # Calculate the center of geometry and principal axis for both helices
        cog1 = helix1.center_of_geometry()
        cog2 = helix2.center_of_geometry()
        axis1 = calculate_principal_axis(helix1)
        axis2 = calculate_principal_axis(helix2)
        
        # Calculate helix-helix distance
        distance = np.linalg.norm(np.cross(cog2 - cog1, axis1)) / np.linalg.norm(axis1)
        
        # Calculate crossing angle Ω
        h = (cog2 - cog1) / np.linalg.norm(cog2 - cog1)
        cross_angle = np.degrees(np.arccos(np.dot(np.cross(axis1, h), np.cross(h, axis2)) / 
                                           (np.linalg.norm(np.cross(axis1, h)) * np.linalg.norm(np.cross(h, axis2)))))
        
        # Collect the results
        results.append((mobile.trajectory.time, distance, cross_angle))

    # Convert results to a numpy array
    results = np.array(results)
    
    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, results, fmt='%10.5f %10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--helix1-start', type=int, required=True,
                        help="Starting residue number of the first helix")
    parser.add_argument('--helix1-end', type=int, required=True,
                        help="Ending residue number of the first helix")
    parser.add_argument('--helix2-start', type=int, required=True,
                        help="Starting residue number of the second helix")
    parser.add_argument('--helix2-end', type=int, required=True,
                        help="Ending residue number of the second helix")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_helix_distance_crossing_angle(**settings)


if __name__ == '__main__':
    main()
