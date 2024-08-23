import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.analysis.align import AlignTraj, AverageStructure  # type: ignore

ANALYSIS_NAME = 'rmsd'


def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time RMSD"
    else:
        header_str = "#time RMSD"

    print(header_str, file=outfile)

    return header_str


def write_rmsd(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
               out: sta.FileRef, align_out: io.TextIOWrapper, ref_psf: Optional[sta.FileRef] = None,
               ref_coor: Optional[sta.FileRef] = None,
               ref_frame_type: str = 'specific', ref_frame_num: int = 1,
               interval: int = 1) -> None:
    """Writes RMSD to `out` file"""

    if ref_psf is None:
        ref_psf = cast(sta.FileRef, psf)
    if ref_coor is None:
        ref_coor = cast(sta.FileRef, traj[0])

    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, ref_coor)
    if ref_frame_type == 'specific':
        ref.transfer_to_memory(start=ref_frame_num-1, stop=ref_frame_num)
    elif ref_frame_type == 'average':
        averaged = AverageStructure(mobile, ref, select=sel).run()
        ref = averaged.results.universe
    else:
        print(f"unknown reference frame type: '{ref_frame_type}'", file=sys.stderr)
        sys.exit(1)

    alignment = AlignTraj(mobile, ref, filename=align_out.name, select=sel).run()

    # TODO: infer start/step/end from settings
    start = .1
    step = .1
    end = (alignment.rmsd.shape[0] + 1) * step
    times = np.arange(start, end, step)
    output = np.stack([times, alignment.rmsd]).T

    with sta.resolve_file(out, 'w') as outfile:
        # TODO: infer fmt from settings
        np.savetxt(outfile, output, fmt='%10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('-rp', '--ref_psf', '--ref_psf_path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('-rc', '--ref_coor', '--ref_coor_path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="Coordinate file to use for reference, if not same as --traj")
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for RMSD calculation")
    parser.add_argument('-rt', '--ref_frame_type', choices=['specific', 'average'],
                        default='specific', metavar='TYPE',
                        help="specific: use a given frame as the reference; "
                        "average: use average structure")
    parser.add_argument('-rn', '--ref_frame_num', type=int, default=1, metavar='N',
                        help="Frame to use for reference coordinates (default: 1). "
                        "Only meaningful if --ref-frame-type is 'specific'")
    parser.add_argument('--align_out', type=argparse.FileType('w'),
                        metavar='FILE', default=None,
                        help="Write aligned trajectory to this path")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_rmsd(**settings)


if __name__ == '__main__':
    main()
