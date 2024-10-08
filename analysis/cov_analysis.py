import argparse
from typing import Optional
import matplotlib.pyplot as plt
import numpy as np

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda  # type: ignore
from MDAnalysis.analysis import align, pca
from sklearn.decomposition import PCA

ANALYSIS_NAME = 'cov_analysis'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#correlation matrix= Atomic x,y,z times atomic x,y,z coordinates"

    print(header_str, file=outfile)

    return header_str


def write_correlation_matrix(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel: str,
                             time_step: float | str, interval: int = 1, 
                             corr_matrix_out: sta.FileRef='', 
                             eigenvalues_out: sta.FileRef='',
                             eigenvectors_out: sta.FileRef='') -> None:
    """Writes correlation matrix to `out` file"""

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    step_num = 1

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
    
    atoms = u.select_atoms(sel) # default 'name CA'
    ref = u.select_atoms(sel)  

    aligner = align.AlignTraj(u, ref, select=sel, match_atoms=True, in_memory=True).run()
    
    # Create a DCD writer object to write the aligned trajectory
    # with mda.Writer("align_traj.dcd", atoms.n_atoms) as W: 
    #     for ts in u.trajectory: 
    #         W.write(atoms) 

    positions = []
    for ts in u.trajectory:
        positions.append(atoms.positions.flatten())   # only CA atoms by deafult
    
    positions = np.array(positions)
    positions = positions.reshape(len(u.trajectory), len(atoms), 3)
    #print(positions.shape)
    mean_positions = positions.mean(axis=0)
    #print(mean_positions.shape)
    centered_positions = positions - mean_positions
    #print(centered_positions.shape)

    centered_positions_flat = centered_positions.reshape(len(u.trajectory), len(atoms) * 3)
    #print("Centered positions (flattened) shape:", centered_positions_flat.shape)


    covariance_matrix = np.cov(centered_positions_flat, rowvar=False)
    #print(covariance_matrix.shape)

    # Compute and plot the correlation matrix
    correlation_matrix = np.corrcoef(centered_positions_flat, rowvar=False)  # normalizes the coefficients


    with sta.resolve_file(corr_matrix_out, 'w') as outfile1:
        header(outfile1)
        np.savetxt(outfile1, correlation_matrix, fmt='%.6f')

    correlation_matrix = np.loadtxt('correlation_matrix.dat', skiprows=1)
    plt.figure(figsize=(10, 8))
    plt.imshow(correlation_matrix, cmap='coolwarm', interpolation='none')
    plt.colorbar(label='Correlation')
    plt.title('Correlation Matrix')
    plt.xlabel('Atomic x,y,z Coordinates')
    plt.ylabel('Atomic x,y,z Coordinates')
    plt.xlim(0, atoms.n_atoms * 3)
    plt.ylim(0, atoms.n_atoms * 3)
    plt.grid()
    plt.savefig('correlation_matrix_heatmap.png', format='png', dpi=300)

    # Perform PCA
    pca = PCA()
    pca.fit(covariance_matrix)

    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_variance_sum = np.cumsum(explained_variance_ratio)
    # Find the number of components that explain at least threshold of the variance
    threshold=0.85
    num_components = np.argmax(cumulative_variance_sum >= threshold) + 1
    #print(f"Number of principal components explaining {threshold * 100}% of the variance: {num_components}")


    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

    # Sort eigenvalues and eigenvectors in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    # Save eigenvalues to file
    with sta.resolve_file(eigenvalues_out, 'w') as outfile2:
        np.savetxt(outfile2, eigenvalues[:num_components], fmt='%.6f', header='Top Eigenvalues')

    # Save eigenvectors to file
    with sta.resolve_file(eigenvectors_out, 'w') as outfile3:
        np.savetxt(outfile3, eigenvectors[:num_components], fmt='%.6f', header='Top Eigenvalues')

    #Print top eigenvalues and eigenvectors
    #print("Top Eigenvalues:")
    #print(eigenvalues[:num_components])  # should be printed as output
    #print("Top Eigenvectors:")
    #print(eigenvectors[:num_components]) # should be printed as output


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
#    parser.add_argument('-c', '--center', action='store_true')
    parser.add_argument('--sel', metavar='selection'),
    parser.add_argument('--corr_matrix_out', metavar='selection'),
    parser.add_argument('--eigenvalues_out', metavar='selection'),
    parser.add_argument('--eigenvectors_out', metavar='selection',
                        help="Membrane headgroup atom selection")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_correlation_matrix(**settings)


if __name__ == '__main__':
    main()
