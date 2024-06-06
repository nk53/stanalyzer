import argparse
import sys
from typing import Optional

from stanalyzer.bin.stanalyzer import FileLike, add_project_args, get_settings, from_settings
import MDAnalysis as mda    # type: ignore

ANALYSIS_NAME = 'system_size'


def header(outfile: Optional[FileLike] = None, include_angles: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    if include_angles:
        header_str = "#time xtla xtlb xtlc alpha beta gamma volume"
    else:
        header_str = "#time xtla xtlb xtlc volume"

    print(header_str, file=outfile)

    return header_str


def write_system_size(psf: FileLike | str, traj: list[FileLike | str],
                      out: FileLike | str, time_step: float | str,
                      interval: int = 1, include_angles: bool = False) -> None:
    """Writes system size to `out` file"""
    n_fields = 8 if include_angles else 5
    output_fmt = ' '.join(["{:.2f}"]*n_fields)

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    sim_time = time_step
    step_num = 1

    if isinstance(out, str):
        out = open(out, 'w')

    with out as outfile:
        header(outfile, include_angles)
        for traj_file in traj:
            u = mda.Universe(psf, traj_file)
            for ts in u.trajectory:
                if step_num % interval:
                    step_num += 1
                    continue

                # get X/Y/Z and possibly also alpha/beta/gamma
                dim_fields = ts.dimensions if include_angles else ts.dimensions[:3]

                # calc volume
                x, y, z = dim_fields[:3]
                volume = x*y*z

                # write output
                output = output_fmt.format(sim_time, *dim_fields, volume)
                print(output, file=outfile)

                # update time
                sim_time += time_step
                step_num += 1


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='stanalyzer system_size')
    add_project_args(parser, 'psf', 'traj')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout,
                        help="File to write results (default: stdout)")
    parser.add_argument('-ts', '--time_step', type=float, default=from_settings,
                        help="Amount of time between frames in trajectory files")
    parser.add_argument('-i', '--interval', type=int, default=from_settings)
    parser.add_argument('-a', '--include_angles', action='store_true')

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(get_settings(ANALYSIS_NAME))

    write_system_size(**settings)


if __name__ == '__main__':
    main()
