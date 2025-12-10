import argparse
from typing import Optional

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
    # https://scikit-learn-extra.readthedocs.io/en/stable/generated/sklearn_extra.cluster.KMedoids.html
    from sklearn_extra.cluster import KMedoids  # type: ignore[import-not-found]
    u = mda.Universe(psf, traj)
    matrix = diffusionmap.DistanceMatrix(u, select=sel).run()
    dist_matrix = matrix.results.dist_matrix
    kmedoids = KMedoids(n_clusters=k, metric='precomputed').fit(dist_matrix)
    print(kmedoids.labels_)  # Labels of each point
    print(kmedoids.medoid_indices_)  # The indices of the medoid rows in X


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    # TODO: selection
    # parser.add_argument('--sel', metavar='selection',
    #                     help="Atom selection for RMSD calculation")
    parser.add_argument('--k', type=float, help='#cluster')
    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    sel = 'name CA'
    k = 2
    run_clustering(settings['psf'], settings['traj'], sel, k)


if __name__ == '__main__':
    main()
