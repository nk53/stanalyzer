import argparse
import typing as t

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from sklearn.decomposition import PCA

import stanalyzer.cli.stanalyzer as sta

ANALYSIS_NAME = 'cov_analysis'


def header(outfile: t.Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#correlation matrix= Atomic x,y,z times atomic x,y,z coordinates"

    print(header_str, file=outfile)

    return header_str


def write_correlation_matrix(psf: sta.FileRef, traj: sta.FileRefList,
                             out: sta.FileRef, sel: str, time_step: float | str,
                             corr_matrix_out: str, eigenvalues_out: str,
                             eigenvectors_out: str, interval: int = 1) -> None:
    """Writes correlation matrix to `out` file"""

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)

    atoms = u.select_atoms(sel)  # default 'name CA'
    ref = u.select_atoms(sel)

    align.AlignTraj(u, ref, select=sel, match_atoms=True, in_memory=True).run()

    ts_positions: list[np.ndarray] = []
    for ts in u.trajectory:
        # only CA atoms by deafult
        ts_positions.append(atoms.positions.flatten())

    positions = np.array(ts_positions)
    positions = positions.reshape(len(u.trajectory), len(atoms), 3)

    mean_positions = positions.mean(axis=0)
    centered_positions = positions - mean_positions
    centered_positions_flat = centered_positions.reshape(
        len(u.trajectory), len(atoms) * 3)

    covariance_matrix = np.cov(centered_positions_flat, rowvar=False)

    # Compute and plot the correlation matrix
    correlation_matrix = np.corrcoef(
        centered_positions_flat, rowvar=False)  # normalizes the coefficients

    with sta.resolve_file(corr_matrix_out, 'w') as outfile1:
        header(outfile1)
        np.savetxt(outfile1, correlation_matrix, fmt='%.6f')

    correlation_matrix = np.loadtxt(corr_matrix_out, skiprows=1)
    plt.figure(figsize=(10, 8))
    plt.imshow(correlation_matrix, cmap='coolwarm', interpolation='none')
    plt.colorbar(label='Correlation')
    plt.title('Correlation Matrix')
    plt.xlabel('Atomic x,y,z Coordinates')
    plt.ylabel('Atomic x,y,z Coordinates')
    plt.grid()
    plt.savefig('correlation_matrix_heatmap.png', format='png', dpi=300)

    # Perform PCA
    pca = PCA()
    pca.fit(covariance_matrix)

    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_variance_sum = np.cumsum(explained_variance_ratio)

    # Find the number of components that explain at least threshold of the variance
    threshold = 0.85
    num_components = np.argmax(cumulative_variance_sum >= threshold) + 1

    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

    # Sort eigenvalues and eigenvectors in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    # Save eigenvalues to file
    with sta.resolve_file(eigenvalues_out, 'w') as outfile2:
        np.savetxt(outfile2, eigenvalues[:num_components],
                   fmt='%.6f', header='Top Eigenvalues')

    # Save eigenvectors to file
    with sta.resolve_file(eigenvectors_out, 'w') as outfile3:
        np.savetxt(outfile3, eigenvectors[:num_components],
                   fmt='%.6f', header='Top Eigenvalues')


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('--sel', metavar='selection', required=True),
    parser.add_argument('--corr-matrix-out', metavar='filename', required=True),
    parser.add_argument('--eigenvalues-out', metavar='filename', required=True),
    parser.add_argument('--eigenvectors-out', metavar='filename', required=True)

    return parser


def main(settings: t.Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_correlation_matrix(**settings)


if __name__ == '__main__':
    main()
