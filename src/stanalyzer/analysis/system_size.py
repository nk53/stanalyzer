import argparse

import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta

ANALYSIS_NAME = 'system_size'


def header(outfile: sta.FileLike | None = None, include_angles: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    if include_angles:
        header_str = "#time xtla xtlb xtlc alpha beta gamma volume"
    else:
        header_str = "#time xtla xtlb xtlc volume"

    print(header_str, file=outfile)

    return header_str


def write_system_size(psf: sta.FileRef, traj: sta.FileRefList,
                      out: sta.FileRef, time_step: float | str,
                      interval: int = 1, include_angles: bool = False) -> None:
    """Writes system size to `out` file"""
    n_fields = 8 if include_angles else 5
    output_fmt = ' '.join(["{:0<#10.5g}"]*n_fields)

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    sim_time = time_step

    with sta.resolve_file(out) as outfile:
        header(outfile, include_angles)
        for traj_file in traj:
            u = mda.Universe(psf, traj_file)
            for step_num, ts in enumerate(u.trajectory, start=1):
                if step_num % interval:
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


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'time_step', 'interval')
    parser.add_argument('-a', '--include-angles', action='store_true')

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_system_size(**settings)


if __name__ == '__main__':
    main()
