import argparse
import sys
from typing import Optional
import io

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.transformations import center_in_box  # type: ignore
from MDAnalysis.core.groups import AtomGroup  # type: ignore

ANALYSIS_NAME = 'position_time'

def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file."""
    if np_formatted:
        header_str = "time input_axis"
    else:
        header_str = "#time input_axis"
    print(header_str, file=outfile)
    return header_str


def write_position_time(psf: sta.FileRef, traj: sta.FileRefList,
                        out: sta.FileRef, sel: str,
                        head_group: str,
                        method: Optional[str] = None,
                        axis: Optional[str] = None) -> None:

    # Load the traj and topology
    u = mda.Universe(psf, traj)

    # Pass the atoms for membrane and the ligand (or selected group for time series calculation) 
    selected_atoms = u.select_atoms(sel)
    selected_head_group = u.select_atoms(head_group)

    # Prepare time series list
    time_series = []

    # Mapping axis (x, y, or z) to the appropriate index for coordinates
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_index = axis_map[axis]  # Use axis parameter from arguments

    # Iterate through each frame in the trajectory
    for ts in u.trajectory:
        # Center the membrane 
        ts = center_in_box(selected_head_group, point=(0, 0, 0))(ts)
        
        # Centroid of the selected atoms (e.g., ligand)
        if method == 'com':
            centroid = selected_atoms.center_of_mass()
        else:
            centroid = selected_atoms.center_of_geometry()

        # Append time and selected axis position of the ligand's centroid to the time series
        time_series.append([u.trajectory.time, centroid[axis_index]])

    # u = mda.Universe(psf, traj)
    # selected_atoms = u.select_atoms(sel)
    # time_series = []

    # axis_map = {'x': 0, 'y': 1, 'z': 2}
    # axis_index = axis_map[axis]  # Use axis parameter from arguments

    # for ts in u.trajectory:
    #     if method == 'com':
    #         centroid = selected_atoms.center_of_mass()
    #     else:
    #         centroid = selected_atoms.center_of_geometry()

    #     time_series.append([u.trajectory.time, centroid[axis_index]])
        
    # print("printing", axis_index, "position")
    # print(time_series)
    
    # Save to output file
    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, time_series, fmt='%10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')

    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection to print position for")
    parser.add_argument('-method', '--method', choices=['com', 'cog'],
                        default='com', metavar='TYPE')
    parser.add_argument('-axis', '--axis', choices=['x', 'y', 'z'],
                        default='z', metavar='TYPE')
    parser.add_argument('--head-group', metavar='selection',
                        help="Membrane headgroup atom selection")
    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_position_time(**settings)


if __name__ == '__main__':
    main()
