import argparse

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda  # type: ignore
from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis as WBA  # type: ignore

ANALYSIS_NAME = 'water_bridge'


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#frame site1 water_to_site1 water_to_site2 site2 (each site and water are in " +\
                 "the form of [acceptor, None] or [donor_hydrogen, donor_heavy])"

    print(header_str, file=outfile)

    return header_str


def print_network(universe, frame, network, out_lines, water_indices):
    if network is None:
        return
    else:
        for node in network:
            if network[node] is None:
                continue
            # Search site1...water_to_site1...water_to_site2...site2
            # node contains site1...water_to_site1
            # subnode contains water_to_site2...site2
            if node[0] in water_indices or node[2] not in water_indices:
                # skip if site1 is water or water_to_site1 is not water
                continue
            for subnode in network[node]:
                if subnode[0] not in water_indices or subnode[2] in water_indices:
                    # skip if water_to_site2 is not water or site2 is water
                    continue
                out_line = str(frame)
                for i in range(4):
                    if i == 0:
                        out_line += ' ['
                    elif i == 2:
                        out_line += '] ['
                    else:
                        out_line += ' '
                    if node[i] is None:
                        out_line += 'None'
                    else:
                        atom = universe.atoms[node[i]]
                        out_line += '_'.join(atom.segid, atom.resname, str(atom.resid), atom.name)
                for i in range(4):
                    if i == 0 or i == 2:
                        out_line += '] ['
                    else:
                        out_line += ' '
                    if subnode[i] is None:
                        out_line += 'None'
                    else:
                        atom = universe.atoms[subnode[i]]
                        out_line += '_'.join(atom.segid, atom.resname, str(atom.resid), atom.name)
                out_line += ']'
                out_lines.append(out_line)
            print_network(universe, frame, network[node], out_lines, water_indices)


def write_water_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                       sel: str = 'protein', sel2: str | None = None, water_sel: str | None = None,
                       d_a_cutoff: float = 3.0, d_h_a_angle_cutoff: float = 150.0,
                       interval: int = 1) -> None:
    """Writes water bridge to `out` file"""

    u = mda.Universe(psf, traj)

    if sel2 is None or sel2.lower() == 'none':
        sel2 = sel
    if water_sel is None or water_sel == 'none':
        water_sel = 'resname TIP3* TP3* TIP4* HOH HO4 WAT SOL'

    water = u.select_atoms(water_sel)
    water_indices = set(water.indices)

    w = WBA(universe=u, selection1=sel, selection2=sel2, water_selection=water_sel,
            order=1,  # consider the pairs connected by one water molecule only
            distance=d_a_cutoff, angle=d_h_a_angle_cutoff)
    w.run(step=interval)

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        out_lines: list[str] = []
        for i in range(len(w.results.network)):
            print_network(u, i*interval, w.results.network[i], out_lines, water_indices)
        for line in out_lines:
            print(line, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection', default='protein',
                        help="Restrict the search of atoms that are bridged by water to only those "
                             "atoms")
    parser.add_argument('--sel2', metavar='selection', default=None,
                        help="Optional. The second group of atoms when searching water bridges "
                             "between two different groups")
    parser.add_argument('--water-sel', metavar='selection', default=None,
                        help="Atom selection for bridging water. If None, then all water molecules "
                             "will be selected")
    parser.add_argument('--d-a-cutoff', type=float, metavar='N', default='3.0',
                        help="Distance cutoff for hydrogen bonds. This cutoff refers to the D-A "
                        "distance.")
    parser.add_argument('--d-h-a-angle-cutoff', type=float, metavar='N', default='150.0',
                        help="D-H-A angle cutoff for hydrogen bonds (degrees).")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_water_bridge(**settings)


if __name__ == '__main__':
    main()
