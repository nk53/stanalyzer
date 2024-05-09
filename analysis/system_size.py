from io import TextIOWrapper
from typing import List, Optional

from stanalyzer.bin import get_settings
import MDAnalysis as mda    # type: ignore

ANALYSIS_NAME = 'system_size'


def header(outfile: Optional[TextIOWrapper] = None, include_angles: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file
    """
    if include_angles:
        header_str = "#time xtla xtlb xtlc alpha beta gamma volume"
    else:
        header_str = "#time xtla xtlb xtlc volume"

    print(header_str, file=outfile)

    return header_str


def write_system_size(psf: str, traj: List[str], out: str, time_step: float, interval: int = 1,
                      include_angles: bool = False) -> None:
    """Writes system size to `out` file
    """
    n_fields = 8 if include_angles else 5
    output_fmt = ' '.join(["{:.2f}"]*n_fields)

    sim_time = time_step
    step_num = 1

    with open(out) as outfile:
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


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = get_settings(ANALYSIS_NAME)
    write_system_size(**settings)


if __name__ == '__main__':
    main()
