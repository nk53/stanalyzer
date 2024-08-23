import argparse
from typing import Optional
import numpy as np

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda  # type: ignore

ANALYSIS_NAME = 'density_z'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#bin density_z"

    print(header_str, file=outfile)

    return header_str


def write_density(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel: str,
                  time_step: float | str, center: bool = False, interval: int = 1,
                  bin_size: float | str = 1.0) -> None:
    """Writes density to `out` file"""
    n_fields = 2
    output_fmt = ' '.join(["{:0<#10.5g}"]*n_fields)

    if isinstance(bin_size, str):  # converting bin_size into float if read as string
        bin_size = float(bin_size.split()[0])   # default value is 1.0

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    step_num = 1

    zc = []  # list of z coordiantes for selected atoms
    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
        head_group = u.select_atoms(sel)
        for ts in u.trajectory:
            if step_num % interval == 0:
                z_cor = u.atoms.positions[:, 2]  # done on all the atoms
                zc.append(head_group.positions[:, 2])  # appended the list for selected atoms
                step_num += 1
                continue

    # computed over all the atoms
    z_min = z_cor.min()
    z_max = z_cor.max()

    # defining grid_widths last value is not counted
    nbins = np.arange(z_min, z_max + bin_size, bin_size)

    density_profile = np.zeros(len(nbins)-1)  # array

    # density is written as a function of edges by default
    density, edges = np.histogram(zc, bins=nbins)

    norm = sum(density)

    bin_centers = [z_min + i*bin_size for i in range(len(nbins)-1)]  # list

    density_profile = density/norm
    density_profile = density_profile.tolist()  # converting array to list

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)

        for bin_center, density in zip(bin_centers, density_profile):
            output = output_fmt.format(bin_center, density)
            print(output, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('-c', '--center', action='store_true')
    parser.add_argument('--sel', metavar='selection'),
    parser.add_argument('--bin_size', metavar='selection', type=float,
                        help="Membrane headgroup atom selection")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_density(**settings)


if __name__ == '__main__':
    main()
