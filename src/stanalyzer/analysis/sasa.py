import argparse
import os

import MDAnalysis as mda  # type: ignore
import freesasa  # type: ignore

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

ANALYSIS_NAME = 'sasa'


def header(outfile: sta.FileLike | None = None, np_formatted=False) -> str:
    """Returns a header string and, optionally writes it to a file."""
    if np_formatted:
        header_str = "frame SASA"
    else:
        header_str = "# frame SASA"

    print(header_str, file=outfile)

    return header_str


def write_sasa(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
               probe_radius: float, algorithm: str, out: sta.FileRef,
               interval: int = 1) -> None:
    """Calculate and write SASA to `out` file."""

    u = mda.Universe(psf, traj)
    atom_group = u.select_atoms(sel)

    if algorithm == 'shrake':
        parameters = freesasa.Parameters(
            {'algorithm': freesasa.ShrakeRupley, 'probe-radius': probe_radius})
    elif algorithm == 'lee':
        parameters = freesasa.Parameters(
            {'algorithm': freesasa.LeeRichards, 'probe-radius': probe_radius})
    else:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. Use 'shrake' or 'lee'.")

    with sta.resolve_file(out, 'w') as f_out:
        f_out.write("# Frame    SASA (Å^2)\n")

        for ts in u.trajectory[::interval]:
            temp_pdb = f"temp_selected_frame_{ts.frame}.pdb"
            atom_group.write(temp_pdb)

            structure = freesasa.Structure(temp_pdb)
            result = freesasa.calc(structure, parameters)

            total_sasa = result.totalArea()
            f_out.write(f"{ts.frame}    {total_sasa:.5f}\n")
            print(f"Frame {ts.frame}: SASA = {total_sasa:.5f} Å^2")

            os.remove(temp_pdb)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for SASA calculation")
    parser.add_argument('--probe-radius', type=p_float,
                        default=1.4, help='Probe radius in Ångstroms')
    parser.add_argument('--algorithm', default='shrake',
                        choices=['shrake', 'lee'], help='Algorithm for SASA calculation')
    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_sasa(**settings)


if __name__ == '__main__':
    main()
