import argparse
import io

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.stanalyzer import writable_outfile
import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis.align import rotation_matrix


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

    if interval < 1:
        raise ValueError("interval must be at least 1")

    if ref_psf is None:
        ref_psf = psf

    # Load mobile and reference universes. The reference uses the first
    # trajectory frame, matching the previous AlignTraj implementation.
    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, traj[0])

    mobile_align = mobile.select_atoms(sel_align)
    reference_align = ref.select_atoms(sel_align)
    rmsf_atoms = mobile.select_atoms(sel_rmsf)

    if len(mobile_align) == 0:
        raise ValueError(f"Alignment selection matched no atoms: {sel_align}")
    if len(rmsf_atoms) == 0:
        raise ValueError(f"RMSF selection matched no atoms: {sel_rmsf}")
    if len(mobile_align) != len(reference_align):
        raise ValueError(
            "Mobile and reference alignment selections must contain "
            "the same number of atoms"
        )

    reference_center = reference_align.positions.mean(axis=0)
    reference_coordinates = (
        reference_align.positions.astype(np.float64, copy=True)
        - reference_center
    )

    mean = np.zeros((len(rmsf_atoms), 3), dtype=np.float64)
    sumsquares = np.zeros_like(mean)
    frame_count = 0

    writer = None
    if align_out is not None:
        writer = mda.Writer(align_out.name, n_atoms=len(mobile.atoms))

    try:
        for ts in mobile.trajectory[::interval]:
            mobile_center = mobile_align.positions.mean(axis=0)
            mobile_coordinates = (
                mobile_align.positions.astype(np.float64, copy=False)
                - mobile_center
            )
            rotation, _ = rotation_matrix(
                mobile_coordinates,
                reference_coordinates,
            )

            aligned_positions = (
                (rmsf_atoms.positions - mobile_center) @ rotation.T
                + reference_center
            )

            frame_count += 1
            delta = aligned_positions - mean
            mean += delta / frame_count
            sumsquares += delta * (aligned_positions - mean)

            if writer is not None:
                mobile.atoms.translate(-mobile_center)
                mobile.atoms.rotate(rotation)
                mobile.atoms.translate(reference_center)
                writer.write(mobile.atoms)
    finally:
        if writer is not None:
            writer.close()

    if frame_count == 0:
        raise ValueError("Trajectory contains no frames")

    rmsf = np.sqrt(sumsquares.sum(axis=1) / frame_count)
    output = np.stack([rmsf_atoms.resids, rmsf]).T

    # Write the results to the output file
    with sta.resolve_file(out, 'w') as outfile:
        np.savetxt(outfile, output, fmt='%10.5f %10.5f', header="residue_indices RMSF")

    if align_out is not None:
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
