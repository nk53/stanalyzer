import argparse
import numpy as np
import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta

ANALYSIS_NAME = 'system_size'


def header(outfile: sta.FileLike | None = None, include_angles: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    # if include_angles:
    #     header_str = "#time xtla xtlb xtlc alpha beta gamma volume"
    # else:
    #     header_str = "#time xtla xtlb xtlc volume"
    header_str = "#time xtla xtlb xtlc alpha beta gamma volume"

    print(header_str, file=outfile)

    return header_str


def write_system_size(psf: sta.FileRef, traj: sta.FileRefList,
                      out: sta.FileRef, time_step: float | str,
                      interval: int = 1, include_angles: bool = False) -> None:
    """Writes system size to `out` file"""
    # n_fields = 8 if include_angles else 5
    n_fields = 8  # always include angles
    output_fmt = ' '.join(["{:0<#10.5g}"]*n_fields)

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    sim_time = time_step

    pi = np.acos(-1.0)    # radian
    rad_to_deg = 180.0/pi  # conversion factor from radian to degree

    with sta.resolve_file(out) as outfile:
        header(outfile, include_angles)
        for traj_file in traj:
            u = mda.Universe(psf, traj_file)
            for step_num, ts in enumerate(u.trajectory, start=1):
                if step_num % interval:
                    sim_time += time_step
                    continue

                # get X/Y/Z and possibly also alpha/beta/gamma
                # dim_fields = ts.dimensions if include_angles else ts.dimensions[:3]
                dim_fields = ts.dimensions

                # calc volume
                x, y, z = dim_fields[:3]
                cos = np.cos(dim_fields[3:]/rad_to_deg)
                cos2 = cos * cos
                cosabc = cos[0] * cos[1] * cos[2]
                # print(cos,cos2,cosabc,1-np.sum(cos2)+2.0*cosabc)
                volume = x*y*z * np.sqrt(1.0 - np.sum(cos2) + 2.0 * cosabc)
                # print(x,y,z,dim_fields[3],dim_fields[4],dim_fields[5],volume)
                # sys.exit(0)

                # write output
                output = output_fmt.format(sim_time, *dim_fields, volume)
                print(output, file=outfile)

                # update time
                sim_time += time_step


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'time_step', 'interval')

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_system_size(**settings)


if __name__ == '__main__':
    main()
