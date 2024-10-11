import argparse
from typing import Optional

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda  # type: ignore
from MDAnalysis.analysis import align, pca  # type: ignore

ANALYSIS_NAME = 'cov_analysis'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#correlation matrix= Atomic x,y,z times atomic x,y,z coordinates"

    print(header_str, file=outfile)

    return header_str


def write_correlation_matrix(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel: str,
                             time_step: float | str, interval: int = 1,
                             sel_atoms_pdb_out: sta.FileRef = '',
                             align_traj_out: sta.FileRef = '',
                             corr_matrix_out: sta.FileRef = '',
                             projected_traj_out: sta.FileRef = '') -> None:
    """Writes correlation matrix, selection.pdb, align_traj.dcd, projected_traj.dcd to `out` file"""

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)

    atoms = u.select_atoms(sel)  # default 'name CA'

    align.AlignTraj(u, u, sel, match_atoms=True, in_memory=True).run()

    # Create a pdb for selection atoms
    with mda.Writer(sel_atoms_pdb_out, atoms.n_atoms) as f:
        u.trajectory[0]
        f.write(atoms)

    # Aligned trajectory for selected atoms to a new DCD file
    with mda.Writer(align_traj_out, atoms.n_atoms) as f:
        for ts in u.trajectory:
            f.write(atoms)

    ts_positions = []
    for ts in u.trajectory:
        ts_positions.append(atoms.positions.flatten())

    positions = np.array(ts_positions)
    positions = positions.reshape(len(u.trajectory), len(atoms), 3)
    mean_positions = positions.mean(axis=0)
    centered_positions = positions - mean_positions

    centered_positions_flat = centered_positions.reshape(
        len(u.trajectory), len(atoms) * 3)

    # Compute correlation matrix
    correlation_matrix = np.corrcoef(centered_positions_flat, rowvar=False)

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
    pc = pca.PCA(u, sel, align=True, mean=None, n_components=None).run()
    ca = u.select_atoms(sel)

    transformed = pc.transform(ca, n_components=3)

    pc1 = pc.p_components[:, 0]
    trans1 = transformed[:, 0]

    projected = np.outer(trans1, pc1) + pc.mean.flatten()
    coordinates = projected.reshape(len(trans1), -1, 3)

    proj1 = mda.Merge(ca)  # new universe for selected atoms
    proj1.load_new(coordinates, order='fac')

    n_frames = transformed.shape[0]

    with mda.Writer(projected_traj_out, atoms.n_atoms) as f:
        for i in range(n_frames):
            coords = coordinates[i]
            proj1.load_new(coords)
            f.write(proj1.atoms)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('--sel', metavar='selection'),
    parser.add_argument('--sel_atoms_pdb_out', metavar='selection'),
    parser.add_argument('--align_traj_out', metavar='selection'),
    parser.add_argument('--corr_matrix_out', metavar='selection'),
    parser.add_argument('--projected_traj_out', metavar='selection',
                        help="Membrane headgroup atom selection")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_correlation_matrix(**settings)


if __name__ == '__main__':
    main()
