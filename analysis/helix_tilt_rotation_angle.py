import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

ANALYSIS_NAME = 'helix_tilt_rotation_angle'

def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time Tilt_Angle Rotation_Angle"
    else:
        header_str = "#time Tilt_Angle Rotation_Angle"

    print(header_str, file=outfile)

    return header_str


def write_tilt_rotation_angle(psf: sta.FileRef, traj: sta.FileRefList,
                              helix_start: int, helix_end: int, out: sta.FileRef, interval: int = 1) -> None:
    """Writes Tilt and Rotation Angles to `out` file"""

    # Read traj and psf
    mobile = mda.Universe(psf, traj)
    
    # Select helix atoms based on residue range (helix_start to helix_end) within the selection `sel`
    helix = mobile.select_atoms(f"resid {helix_start}-{helix_end} and name CA")


    # Calculate the principal axis
    def calculate_principal_axis(atoms):
        positions = atoms.positions - atoms.center_of_geometry()
        inertia_tensor = np.dot(positions.T, positions)
        eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
        principal_axis = eigenvectors[:, np.argmin(eigenvalues)]
        return principal_axis / np.linalg.norm(principal_axis)
    
    # Calculate the tilt angle
    def calculate_tilt_angle(axis):
        z_axis = np.array([0, 0, 1])
        tilt_angle = np.degrees(np.arccos(np.dot(axis, z_axis)))
        return tilt_angle

    # Calculate the rotation angle ρ
    def calculate_rotation_angle(atoms, axis):
        positions = atoms.positions - atoms.center_of_geometry()
        reference_vector = np.array([1, 0, 0])  # Example reference vector in the x-direction
        angles = np.arctan2(np.dot(positions, np.cross(reference_vector, axis)),
                            np.dot(positions, reference_vector))
        return np.degrees(np.mean(angles))

    results = []
    for ts in mobile.trajectory[::interval]:
        axis = calculate_principal_axis(helix)
        tilt_angle = calculate_tilt_angle(axis)
        rotation_angle = calculate_rotation_angle(helix, axis)
        results.append((mobile.trajectory.time, tilt_angle, rotation_angle))

    results = np.array(results)
    
    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, results, fmt='%10.5f %10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--helix_start', type=int, required=True,
                        help="Starting residue number of the helix")
    parser.add_argument('--helix_end', type=int, required=True,
                        help="Ending residue number of the helix")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_tilt_rotation_angle(**settings)


if __name__ == '__main__':
    main()