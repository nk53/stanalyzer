import argparse
from typing import Optional

import kmedoids
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap

import stanalyzer.cli.stanalyzer as sta

ANALYSIS_NAME = 'clustering_kmedoid'


def header(outfile: sta.FileLike | None = None,
           np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time cluster"
    else:
        header_str = "#time cluster"

    print(header_str, file=outfile)

    return header_str


def run_clustering(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   k: int) -> None:
    u = mda.Universe(psf, traj)
    matrix = diffusionmap.DistanceMatrix(u, select=sel).run()
    dist_matrix = matrix.results.dist_matrix
    # fixed random_state for byte-identical golden comparisons
    result = kmedoids.fasterpam(dist_matrix, k, random_state=0)

    # write_to_outfile resolves relpaths to <output_path>/<analysis_name>/
    sout = ''.join(f'{i:8d} {c:8d}\n' for i, c in enumerate(result.labels))
    sta.write_to_outfile('cluster.dat', sout)

    frames = set(result.medoids.tolist())
    protein = u.select_atoms("protein")
    with mda.Writer(sta.writable_outfile('cluster_representative.pdb').name,
                    multiframe=True) as pdb:
        for ts in u.trajectory:
            if ts.frame in frames:
                pdb.write(protein)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    # TODO: selection
    # parser.add_argument('--sel', metavar='selection',
    #                     help="Atom selection for RMSD calculation")
    parser.add_argument('--k', type=int, help='#cluster')
    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    sel = 'name CA'
    k = settings.get('k') or 2
    run_clustering(settings['psf'], settings['traj'], sel, k)


if __name__ == '__main__':
    main()
