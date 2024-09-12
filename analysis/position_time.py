import argparse
import sys
from typing import Optional
import io

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

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
                        method: Optional[str] = None,
                        axis: Optional[str] = None) -> None:

    u = mda.Universe(psf, traj)
    selected_atoms = u.select_atoms(sel)
    time_series = []

    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_index = axis_map[axis]  # Use axis parameter from arguments

    for ts in u.trajectory:
        if method == 'com':
            centroid = selected_atoms.center_of_mass()
        else:
            centroid = selected_atoms.center_of_geometry()

        time_series.append([u.trajectory.time, centroid[axis_index]])
        
    print("printing", axis_index, "position")
    print(time_series)
    
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
    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_position_time(**settings)


if __name__ == '__main__':
    main()