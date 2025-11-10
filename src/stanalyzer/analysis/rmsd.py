import argparse
import sys
import io
from typing import Optional, TypeAlias, cast

import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj, AverageStructure
import numpy as np

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_int

ANALYSIS_NAME = 'rmsd'

OptFileRef: TypeAlias = sta.LazyFile | io.TextIOWrapper | str | None


def header(outfile: Optional[sta.FileLike] = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time RMSD"
    else:
        header_str = "#time RMSD"

    print(header_str, file=outfile)

    return header_str


def write_rmsd(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
               out: sta.FileRef, align_out: OptFileRef,
               ref_psf: Optional[sta.FileRef] = None,
               ref_coor: Optional[sta.FileRef] = None,
               ref_frame_type: str = 'specific', ref_frame_num: int = 1,
               time_step: float | str = .1, interval: int = 1) -> None:
    """Writes RMSD to `out` file"""

    align_filename: str | None

    if ref_psf is None:
        ref_psf = cast(sta.FileRef, psf)
    if ref_coor is None:
        ref_coor = cast(sta.FileRef, traj[0])
    if align_out is None or not align_out:
        align_filename = None
    else:
        align_filename = cast(io.TextIOWrapper, sta.resolve_file(align_out)).name

    time_step = cast(float, sta.resolve_ts(time_step))

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

    alignment = AlignTraj(mobile, ref, filename=align_filename, select=sel).run()

    # TODO: settings for start/end?
    start = time_step
    end = (alignment.rmsd.shape[0] + 1) * time_step
    n_step = mobile.trajectory.n_frames
    times = np.linspace(start, end, n_step)
    output = np.stack([times, alignment.rmsd]).T

    # TODO: avoid calculating RMSD for skipped frames?
    times = times[::interval]
    output = output[::interval]

    with sta.resolve_file(out, 'w') as outfile:
        # TODO: settings for output format?
        np.savetxt(outfile, output, fmt='%10.5f %10.5f', header=header(np_formatted=True))


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('-rp', '--ref-psf', '--ref-psf-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('-rc', '--ref-coor', '--ref-coor-path', type=sta.ExistingFile,
                        metavar='FILE',
                        help="Coordinate file to use for reference, if not same as --traj")
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for RMSD calculation")
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

    write_rmsd(**settings)


if __name__ == '__main__':
    main()
