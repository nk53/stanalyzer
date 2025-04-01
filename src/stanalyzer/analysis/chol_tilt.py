import argparse
import sys
from typing import Optional

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda
from MDAnalysis.transformations import center_in_box
from MDAnalysis.core.groups import AtomGroup
import numpy
import math

ANALYSIS_NAME = 'sterol tilt angle to Z axis'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#time resid degree"

    print(header_str, file=outfile)

    return header_str


def write_chol_degree(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel: str,
                      time_step: float | str, center: bool = False, interval: int = 1) -> None:
    """Writes sterol tilt angle to `out` file"""
    n_fields = 3
    output_fmt = ' '.join(["{:0<#10.5g}"]*n_fields)

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    sim_time = time_step
    step_num = 1

    # sele Sterols # segname MEMB and resname ( CHOL / ERGO / CAMP)
    chol_sel = f"({sel})"

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        for traj_file in traj:
            u = mda.Universe(psf, traj_file)
            chol_atoms = u.select_atoms(sel)
            for ts in u.trajectory:
                if step_num % interval:
                    step_num += 1
                    continue

                # check for alignment/selection problems
                if not chol_atoms.residues:
                    print(f"Error: '{chol_sel}' did not match any atoms",
                          file=sys.stderr)
                    sys.exit(1)

                # membrane centering if selected by user
                if center:
                    origin = 0, 0, 0

                    # picking a single atom ensures pbc boundary doesn't interfere
                    ts = center_in_box(
                        AtomGroup([chol_atoms[0]]), point=origin)(ts)
                    # first centering tends to be a bit off
                    ts = center_in_box(chol_atoms, point=origin)(ts)

                for i in chol_atoms.residues:  # loop through each cholesterol residue
                    C3 = i.atoms.select_atoms('name C3')[0]
                    C17 = i.atoms.select_atoms('name C17')[0]
                    # Sterol tilt angle calculation based on C3 and C17 atoms
                    vector = C3.position - C17.position
                    mag = numpy.linalg.norm(vector)
                    # simple tilt angle calculation to Z axis
                    degree = math.degrees(math.acos(vector[2] / mag))

                    if degree > 90:
                        degree = abs(degree - 180)

                    # write output
                    output = output_fmt.format(sim_time, i.resid, degree)
                    # time resid degree
                    # 1.2 32 57.3
                    # 1.2 36 55.6
                    # 1.2 45 51.7
                    # 1.3 32 57.5
                    # 1.3 36 53.4
                    # 1.3 45 55.9
                    # ... the output should look like this
                    print(output, file=outfile)

                # update time
                sim_time += time_step
                step_num += 1


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('-c', '--center', action='store_true')
    parser.add_argument('--sel', metavar='selection',
                        help="Sterol lipid selection. E.g.: segid MEMB and resname CHL1")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_chol_degree(**settings)


if __name__ == '__main__':
    main()
