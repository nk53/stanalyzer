import argparse
import warnings

import MDAnalysis as mda  # type: ignore
from MDAnalysis.analysis.rdf import InterRDF  # type: ignore
import numpy as np

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

warnings.simplefilter("ignore", category=np.VisibleDeprecationWarning)

ANALYSIS_NAME = 'Radial distribution function'


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#Bin    Frequency"

    print(header_str, file=outfile)

    return header_str


def write_rdf(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel1: str,
              sel2: str, bin_size: float):
    """Writes radial distribution fuction to 'out' file"""

    u = mda.Universe(psf, traj)

    # Determine RDF range and nbins
    dimensions = []
    for ts in u.trajectory:
        dimensions.append(ts.dimensions)
    box_l = np.array(list(zip(*dimensions))[:3]).min()
    rdf_range = box_l / 2.0
    nbins = int(rdf_range / bin_size)
    rdf_range = nbins * bin_size

    sel_atoms1 = u.select_atoms(sel1)
    sel_atom2 = u.select_atoms(sel2)

    rdf = InterRDF(sel_atoms1, sel_atom2, range=(0, rdf_range), nbins=nbins)
    rdf.run()

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)

        for i, data in enumerate(rdf.rdf):
            print(f'{rdf.bins[i]:.2f} {data:.4f}', file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    parser.add_argument('-sel1', metavar='select', default=None,
                        help="Reference selection")
    parser.add_argument('-sel2', metavar='select', default=None,
                        help="Target selection")
    parser.add_argument('-bin_size', metavar='float',
                        default=None, type=p_float, help="Bin size")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_rdf(**settings)


if __name__ == '__main__':
    main()
