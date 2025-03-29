import argparse

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import hole2

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_int

ANALYSIS_NAME = 'Pore Profile'


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#pore radius"

    print(header_str, file=outfile)

    return header_str


def write_pore_radius(psf: sta.FileRef, traj: sta.FileRefList, hist_out: sta.FileRef,
                      midpoints_out: sta.FileRef, means_out: sta.FileRef,
                      sel: str, bins: int = 100, interval: int = 1,
                      hole_path: str = 'hole') -> None:
    """Writes pore radius to out files"""

    traj = [sta.resolve_file(t) for t in traj]
    u = mda.Universe(psf, traj)

    breakpoint()

    with hole2.HoleAnalysis(u, select=sel, cpoint='center_of_geometry',
                            executable=hole_path) as ha:
        ha.run()
        means, edges = ha.histogram_radii(bins=bins, range=None, aggregator=np.mean)

    midpoints = 0.5*(edges[1:]+edges[:-1])

    midpoints = midpoints[::interval]
    means = midpoints[::interval]

    np.savetxt(midpoints_out, midpoints)
    np.savetxt(means_out, means)
    plt.plot(midpoints, means)
    plt.ylabel(r"Mean HOLE radius $R$ ($\AA$)")
    plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
    plt.savefig(hist_out)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval')
    sta.add_exec_args(parser, 'hole')
    parser.add_argument('-s', '--sel', default='protein',
                        help="Atom selection containing pore. Default: protein.")
    parser.add_argument('-b', '--bins', default=100, type=p_int,
                        help="Number of histogram bins in output. Default: 100.")
    parser.add_argument('-ho', '--hist-out', type=argparse.FileType('w'), default="hist.png",
                        help="File to write histogram image. Default: hist.png")
    parser.add_argument('-bo', '--midpoints-out', '--bins-out', type=argparse.FileType('w'),
                        default="midpoints.dat",
                        help="File to write bin midpoints. Default: midpoints.dat")
    parser.add_argument('-mo', '--means-out', type=argparse.FileType('w'), default="means.dat",
                        help="File to write means image. Default: means.dat")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_pore_radius(**settings)


if __name__ == '__main__':
    main()
