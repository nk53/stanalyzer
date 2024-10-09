import argparse
import sys
import io
from typing import Optional, cast

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
from MDAnalysis.analysis import diffusionmap
import numpy as np
from MDAnalysis.analysis.align import AlignTraj
import scipy.cluster.hierarchy as sch

ANALYSIS_NAME = 'clustering_hca'

def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time cluster"
    else:
        header_str = "#time cluster"

    print(header_str, file=outfile)

    return header_str

def run_clustering(psf: sta.FileRef, traj: sta.FileRefList, sel: str, cutoff: float, criterion: str, method: str, write_file: bool, write_file_represent: bool) -> None:
    # https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/pairwise_rmsd.html
    #from MDAnalysisData import datasets
    #from MDAnalysis.tests.datafiles import PSF, DCD
    #u = mda.Universe(PSF, DCD)
    u = mda.Universe(psf, traj)
    aligner = AlignTraj(u, u, select=sel, in_memory=True).run()
    matrix = diffusionmap.DistanceMatrix(u, select=sel).run()
    dist_matrix = matrix.results.dist_matrix

    # datastructure
    clus = None
    clust2index = {}
    represent = {}

    # hierarchical clustering
    n = dist_matrix.shape[0]
    idx = np.triu_indices(n, 1)
    Y = sch.linkage(dist_matrix[idx], method=method) # linkage matrix
    clus = sch.fcluster(Y, cutoff, criterion=criterion)
    num_clus = max(clus)
    for i in range(1, num_clus+1):
        dmat = dist_matrix[clus==i, :]
        dmat = dmat[:, clus==i]
        m = dmat.shape[0]
        index = np.where(clus==i)[0]
        target = list(zip(index, dmat.mean(axis=1)))
        target.sort(key=lambda X: X[1])
        clust2index[i] = index
        represent[i] = target.pop(0)[0]

    # write output
    if write_file:
        f = open('cluster.dat', 'w')
        for i,c in enumerate(clus):
            f.write("%8d %8d\n" % (i,c)) # index, cluster
        f.close()
    if write_file_represent:
        frames = represent.values()
        protein = u.select_atoms("protein")
        with mda.Writer("cluster_representative.pdb", multiframe=True) as pdb:
            for ts in u.trajectory:
                if ts.frame in frames:
                    pdb.write(protein)
    

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    # TODO: selection
    #parser.add_argument('--sel', metavar='selection',
    #                    help="Atom selection for RMSD calculation")
    parser.add_argument('--cutOff',type=float,
                        help='cut-off for clustering')
    return parser

def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    sel = 'name CA'
    cutoff = 2.0
    criterion = 'distance'
    method = 'average'
    write_file = True
    write_file_represent = True
    run_clustering(settings['psf'], settings['traj'], sel, cutoff, criterion, method, write_file, write_file_represent)

if __name__ == '__main__':
    main()
