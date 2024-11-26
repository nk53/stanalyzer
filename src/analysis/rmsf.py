import argparse
from typing import Optional, cast
import io

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.analysis.align import AlignTraj  # type: ignore
from MDAnalysis.analysis.rms import RMSF  # type: ignore

import stanalyzer.cli.stanalyzer as sta


ANALYSIS_NAME = 'rmsf'


def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time RMSF"
    else:
        header_str = "#time RMSF"

    print(header_str, file=outfile)

    return header_str

def write_rmsf(psf: sta.FileRef, traj: sta.FileRefList, sel_align: str, sel_rmsf: str,
                             out: sta.FileRef, align_out: io.TextIOWrapper, 
                             ref_psf: Optional[sta.FileRef] = None,
                             interval: int = 1) -> None:
    """Writes RMSF (Root Mean Square Fluctuation) to `out` file."""

    if ref_psf is None:
        ref_psf = psf

    # Load mobile and reference universes
    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, traj[0])

    # Align the mobile trajectory to the reference based on the selection for alignment
    alignment = AlignTraj(mobile, ref, filename=align_out.name, select=sel_align).run()

    # Load the aligned trajectory from the saved file
    aligned_mobile = mda.Universe(psf, align_out.name)


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
    np.savetxt(out, output, fmt='%10.5f %10.5f', header="residue_indices RMSF")
    print(f"RMSF results saved to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')

    # Add two separate selection arguments, one for alignment and one for RMSF calculation
    parser.add_argument('--sel_align', metavar='selection_align',
                        help="Atom selection for trajectory alignment")
    parser.add_argument('--sel_rmsf', metavar='selection_rmsf',
                        help="Atom selection for RMSF calculation")
    parser.add_argument('--align_out', type=argparse.FileType('w'),
                        metavar='FILE', default=None,
                        help="Write aligned trajectory to this path")

    return parser



def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        parser = get_parser()
        args = parser.parse_args()
        settings = vars(args)

    write_rmsf(**settings)


if __name__ == '__main__':
    main()
