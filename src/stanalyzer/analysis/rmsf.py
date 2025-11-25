import argparse
import io

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.stanalyzer import writable_outfile
import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis.align import AlignTraj
from MDAnalysis.analysis.rms import RMSF


ANALYSIS_NAME = 'rmsf'


def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time RMSF"
    else:
        header_str = "#time RMSF"

    print(header_str, file=outfile)

    return header_str


def write_rmsf(psf: sta.FileRef, traj: sta.FileRefList, sel_align: str, sel_rmsf: str,
               out: sta.FileRef, align_out: io.TextIOWrapper | None = None,
               ref_psf: sta.FileRef | None = None,
               interval: int = 1) -> None:
    """Writes RMSF (Root Mean Square Fluctuation) to `out` file."""

    if ref_psf is None:
        ref_psf = psf

    # Load mobile and reference universes
    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, traj[0])

    align_file = align_out.name if align_out else None

    # Align the mobile trajectory to the reference based on the selection for alignment
    AlignTraj(mobile, ref, filename=align_file, select=sel_align).run()

    # Load the aligned trajectory from the saved file
    aligned_mobile = mobile if align_out is None else mda.Universe(psf, align_file)

    # Calculate RMSF using the aligned trajectory and the selection for RMSF calculation
    # rmsf_analysis = RMSF(mobile.select_atoms(sel_rmsf)).run()
    atoms = aligned_mobile.select_atoms(sel_rmsf)
    rmsf_analysis = RMSF(aligned_mobile.select_atoms(sel_rmsf)).run()
    # atom_index = np.arange(len(rmsf_analysis.rmsf))

    # Combine atom index and RMSF values
    # output = np.stack([atom_index, rmsf_analysis.rmsf]).T
    residue_indices = atoms.resids
    output = np.stack([residue_indices, rmsf_analysis.rmsf]).T

    # Write the results to the output file
    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, output, fmt='%10.5f %10.5f', header="residue_indices RMSF")
    print(f"RMSF results saved to {outfile.name}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')

    # Add two separate selection arguments, one for alignment and one for RMSF calculation
    parser.add_argument('-rp', '--ref-psf', '--ref-psf-path', type=sta.InputFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('--sel-align', metavar='selection_align',
                        help="Atom selection for trajectory alignment")
    parser.add_argument('--sel-rmsf', metavar='selection_rmsf',
                        help="Atom selection for RMSF calculation")
    parser.add_argument('--align-out', type=writable_outfile,
                        metavar='FILE', default=None,
                        help="Write aligned trajectory to this path")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        parser = get_parser()
        args = parser.parse_args()
        settings = vars(args)

    write_rmsf(**settings)


if __name__ == '__main__':
    main()
