import argparse
import sys
from typing import Optional

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda  # type: ignore
from MDAnalysis.transformations import center_in_box  # type: ignore
from MDAnalysis.core.groups import AtomGroup  # type: ignore

ANALYSIS_NAME = 'thickness'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#time thickness"

    print(header_str, file=outfile)

    return header_str


def write_thickness(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef, sel: str,
                    time_step: float | str, center: bool = False, interval: int = 1) -> None:
    """Writes membrane thickness to `out` file"""
    n_fields = 2
    output_fmt = ' '.join(["{:0<#10.5g}"]*n_fields)

    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    sim_time = time_step
    step_num = 1

    up_sel = f"({sel}) and prop z >= 0"
    dn_sel = f"({sel}) and prop z < 0"

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        for traj_file in traj:
            u = mda.Universe(psf, traj_file)
            head_group = u.select_atoms(sel)
            for ts in u.trajectory:
                if step_num % interval:
                    step_num += 1
                    continue

                # membrane centering if selected by user
                if center:
                    origin = 0, 0, 0

                    # picking a single atom ensures pbc boundary doesn't interfere
                    ts = center_in_box(AtomGroup([head_group[0]]), point=origin)(ts)
                    # first centering tends to be a bit off
                    ts = center_in_box(head_group, point=origin)(ts)

                # calculate membrane thickness
                up_atoms = head_group.select_atoms(up_sel)
                dn_atoms = head_group.select_atoms(dn_sel)

                # check for alignment/selection problems
                if not up_atoms:
                    print(f"Error: '{up_sel}' did not match any atoms",
                          file=sys.stderr)
                if not dn_atoms:
                    print(f"Error: '{dn_sel}' did not match any atoms",
                          file=sys.stderr)
                if not (up_atoms and dn_atoms):
                    sys.exit(1)

                _, _, Zup = up_atoms.center_of_geometry()
                _, _, Zdn = dn_atoms.center_of_geometry()

                thickness = Zup - Zdn

                # write output
                output = output_fmt.format(sim_time, thickness)
                print(output, file=outfile)

                # update time
                sim_time += time_step
                step_num += 1


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('-c', '--center', action='store_true')
    parser.add_argument('--sel', metavar='selection',
                        help="Membrane headgroup atom selection")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_thickness(**settings)


if __name__ == '__main__':
    main()
