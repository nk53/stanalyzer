import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.analysis.align import AlignTraj, AverageStructure  # type: ignore

ANALYSIS_NAME = 'radius_of_gyration'


def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time radius_of_gyration"
    else:
        header_str = "#time radius_of_gyration"

    print(header_str, file=outfile)

    return header_str


def write_radius_of_gyration(psf: sta.FileRef, traj: sta.FileRefList, sel_align: str, sel_rg: str,
                             out: sta.FileRef, align_out: io.TextIOWrapper,
                             ref_psf: Optional[sta.FileRef] = None,
                             ref_coor: Optional[sta.FileRef] = None,
                             ref_frame_type: str = 'specific', ref_frame_num: int = 1,
                             interval: int = 1) -> None:
    """Writes the radius of gyration to `out` file"""

    if ref_psf is None:
        ref_psf = cast(sta.FileRef, psf)
    if ref_coor is None:
        ref_coor = cast(sta.FileRef, traj[0])

    # Load mobile and reference universes
    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, ref_coor)

    # Select reference frame
    if ref_frame_type == 'specific':
        ref.transfer_to_memory(start=ref_frame_num-1, stop=ref_frame_num)
    elif ref_frame_type == 'average':
        averaged = AverageStructure(mobile, ref, select=sel_align).run()
        ref = averaged.results.universe
    else:
        print(f"unknown reference frame type: '{ref_frame_type}'", file=sys.stderr)
        sys.exit(1)

    # Align the mobile trajectory to the reference and save the aligned trajectory
    AlignTraj(mobile, ref, filename=align_out.name, select=sel_align).run()

    # Load the aligned trajectory from the saved file
    aligned_mobile = mda.Universe(psf, align_out.name)

    # Calculate radius of gyration for the aligned trajectory
    protein = aligned_mobile.select_atoms(sel_rg)

    Rgyr_list = []
    for ts in aligned_mobile.trajectory[::interval]:  # Iterate over the aligned trajectory
        Rgyr_list.append((aligned_mobile.trajectory.time, protein.radius_of_gyration()))

    Rgyr = np.array(Rgyr_list)

    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, Rgyr, fmt='%10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('-rp', '--ref-psf', '--ref-psf-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('-rc', '--ref-coor', '--ref-coor-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="Coordinate file to use for reference, if not same as --traj")
    parser.add_argument('--sel-align', metavar='selection',
                        help="Atom selection for trajectory alignment")
    parser.add_argument('--sel-rg', metavar='selection',
                        help="Atom selection for radius of gyration calculation")
    parser.add_argument('-rt', '--ref-frame-type', choices=['specific', 'average'],
                        default='specific', metavar='TYPE',
                        help="specific: use a given frame as the reference; "
                             "average: use average structure")
    parser.add_argument('-rn', '--ref-frame-num', type=p_int, default=1, metavar='N',
                        help="Frame to use for reference coordinates (default: 1). "
                             "Only meaningful if --ref-frame-type is 'specific'")
    parser.add_argument('--align-out', type=argparse.FileType('w'),
                        metavar='FILE', default=None,
                        help="Write aligned trajectory to this path")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_radius_of_gyration(**settings)


if __name__ == '__main__':
    main()
