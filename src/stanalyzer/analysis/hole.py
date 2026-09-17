import argparse

import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import hole2

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.stanalyzer import LazyFile, writable_outfile
from stanalyzer.cli.validators import exec_name, p_int

ANALYSIS_NAME = 'Pore Profile'

# HOLE default seeds from wall-clock time
RANDOM_SEED = 42


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#pore radius"

    print(header_str, file=outfile)

    return header_str


def write_pore_radius(psf: sta.FileRef, traj: sta.FileRefList, hist_out: sta.FileRef,
                      midpoints_out: sta.FileRef, means_out: sta.FileRef,
                      sel: str, bins: int = 100, interval: int = 1,
                      hole_path: str | None = 'hole') -> None:
    """Writes pore radius to out files"""

    if hole_path is None:
        # Resolve None → path; raises ValueError if not on PATH
        hole_path = exec_name('hole')
    # mdahole2's HoleAnalysis joins universe/trajectory filenames as str
    traj = [str(t) for t in traj]
    u = mda.Universe(str(psf), traj)

    # CVECT=(0,0,1) prevents HOLE from locking onto a non-pore path
    with hole2.HoleAnalysis(u, select=sel, cpoint='center_of_geometry',
                            cvect=(0.0, 0.0, 1.0), executable=hole_path) as ha:
        # interval strides frames; seed HOLE's Monte-Carlo for reproducibility
        ha.run(step=interval, random_seed=RANDOM_SEED)

        profiles = ha.results.profiles

    with sta.resolve_file(midpoints_out, 'w') as outfile:
        for i, frame in enumerate(sorted(profiles)):
            if i:
                outfile.write('\n')
            profile = profiles[frame]
            np.savetxt(outfile, np.column_stack((profile['rxn_coord'],
                                                 profile['radius'])))

    coords = np.concatenate([profiles[f]['rxn_coord'] for f in profiles])
    radii = np.concatenate([profiles[f]['radius'] for f in profiles])
    edges = np.linspace(coords.min(), coords.max(), bins + 1)
    midpoints = 0.5 * (edges[1:] + edges[:-1])
    counts, _ = np.histogram(coords, bins=edges)
    sums, _ = np.histogram(coords, bins=edges, weights=radii)
    means = np.divide(sums, counts, out=np.zeros(bins), where=counts > 0)
    with sta.resolve_file(means_out, 'w') as outfile:
        np.savetxt(outfile, means)

    plt.plot(midpoints, means)
    plt.ylabel(r"Mean HOLE radius $R$ ($\AA$)")
    plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
    # PNG is binary; pass raw path to matplotlib
    if isinstance(hist_out, LazyFile):
        hist_out = hist_out.name
    plt.savefig(hist_out)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval')
    sta.add_exec_args(parser, 'hole')
    parser.add_argument('-s', '--sel', default='protein',
                        help="Atom selection containing pore. Default: protein.")
    parser.add_argument('-b', '--bins', default=100, type=p_int,
                        help="Number of histogram bins in output. Default: 100.")
    parser.add_argument('-ho', '--hist-out', type=writable_outfile, default="hist.png",
                        help="File to write histogram image. Default: hist.png")
    parser.add_argument('-bo', '--midpoints-out', '--bins-out', type=writable_outfile,
                        default="midpoints.dat",
                        help="File to write bin midpoints. Default: midpoints.dat")
    parser.add_argument('-mo', '--means-out', type=writable_outfile, default="means.dat",
                        help="File to write means image. Default: means.dat")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_pore_radius(**settings)


if __name__ == '__main__':
    main()
