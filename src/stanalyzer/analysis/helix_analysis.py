#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda
import matplotlib.pyplot as plt
from stanalyzer.cli.stanalyzer import writable_outfile
from MDAnalysis.analysis.align import AlignTraj, AverageStructure
from MDAnalysis.analysis import helix_analysis as hel

ANALYSIS_NAME = 'helix_analysis'


def header(outfile: Optional[sta.FileLike] = None,
           np_formatted: bool = False) -> str:
    """Returns a header string and optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time average_twist average_tilt"
    else:
        header_str = "#time average_twist average_tilt"

    print(header_str, file=outfile)

    return header_str


def analyze_helix(aligned: mda.Universe, selection: str,
                  output_file_handle: io.TextIOBase) -> None:

    # Run helix analysis
    h = hel.HELANAL(aligned, select=selection, ref_axis=[0, 0, 1]).run()

    # Extract results
    global_axes = h.results.summary['global_axis']
    global_tilts = h.results.summary['global_tilts']
    all_bends = h.results.all_bends

    # Write results to the file handle
    output_file_handle.write("Global Axes:\n")
    for key, val in global_axes.items():
        output_file_handle.write(f"{key}: {val}\n")

    output_file_handle.write("\nGlobal Tilts:\n")
    for key, val in global_tilts.items():
        output_file_handle.write(f"{key}: {val:.3f}\n")

    output_file_handle.write("\nAll Bends (shape: {}):\n".format(all_bends.shape))
    for i, frame in enumerate(all_bends[:, :, 0]):
        output_file_handle.write(f"Frame {i}: {frame}\n")

    # Optionally plot average local twists per frame
    plt.plot(h.results.local_twists.mean(axis=1))
    plt.xlabel('Frame')
    plt.ylabel('Average twist (degrees)')
    plt.show()


def write_helix_analysis(
        psf: sta.FileRef, traj: sta.FileRefList, sel_align: str,
        sel_helix: str, out: sta.FileRef,
        align_out: io.TextIOWrapper | None = None,
        ref_psf: Optional[sta.FileRef] = None,
        ref_coor: Optional[sta.FileRef] = None,
        ref_frame_type: str = 'specific', ref_frame_num: int = 1,
        interval: int = 1) -> None:
    """Writes the helix analysis results to `out` file"""

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

    align_file = align_out.name if align_out else None

    # Align the mobile trajectory to the reference and save the aligned trajectory
    AlignTraj(mobile, ref, filename=align_file, select=sel_align).run()

    # Load the aligned trajectory from the saved file
    aligned_mobile = mda.Universe(psf, align_file)  # noqa: F841

    # Perform helix analysis and write results
    with sta.resolve_file(out, 'w') as outfile:
        analyze_helix(
            aligned=aligned_mobile,
            selection=sel_helix,
            output_file_handle=outfile  # Pass file handle directly
        )


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
    parser.add_argument('--sel-helix', metavar='selection',
                        help="Atom selection for helix analysis")
    parser.add_argument('-rt', '--ref-frame-type', choices=['specific', 'average'],
                        default='specific', metavar='TYPE',
                        help="specific: use a given frame as the reference; "
                             "average: use average structure")
    parser.add_argument('-rn', '--ref-frame-num', type=int, default=1, metavar='N',
                        help="Frame to use for reference coordinates (default: 1). "
                             "Only meaningful if --ref-frame-type is 'specific'")
    parser.add_argument('--align-out', type=writable_outfile, metavar='FILE',
                        default=None, help="Write aligned trajectory to this path")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_helix_analysis(**settings)


if __name__ == '__main__':
    main()
