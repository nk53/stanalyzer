import argparse
from typing import Optional, cast

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.analysis.align import AlignTraj  # type: ignore
from MDAnalysis.analysis.rms import RMSF  # type: ignore


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


def write_rmsf(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
               out: sta.FileRef, ref_psf: Optional[sta.FileRef] = None,
               interval: int = 1) -> None:
    """Writes RMSF to `out` file"""

    if ref_psf is None:
        ref_psf = cast(sta.FileRef, psf)

    mobile = mda.Universe(psf, traj)
    ref = mda.Universe(ref_psf, traj[0])

    alignment = AlignTraj(mobile, ref, select=sel, in_memory=True).run()

    # Calculate RMSF using the aligned trajectory
    rmsf_analysis = RMSF(mobile.select_atoms(sel)).run()
    atom_index = np.arange(len(rmsf_analysis.rmsf))

    # Combine atom index and RMSF values
    output = np.stack([atom_index, rmsf_analysis.rmsf]).T

    # Write the results to the output file
    np.savetxt(out, output, fmt='%10.5f %10.5f', header="Atom_Index RMSF")
    print(f"RMSF results saved to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for RMSF calculation")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_rmsf(**settings)


if __name__ == '__main__':
    main()
