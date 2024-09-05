import argparse
import sys
from typing import Optional

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

ANALYSIS_NAME = 'salt_bridge'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#frame donor hydrogen acceptor distance angle"

    print(header_str, file=outfile)

    return header_str


#def write_salt_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
#                      sel: str = 'all', hydrogens_sel: str = None, acceptors_sel: str = None, donors_sel: str = None,
#                      d_a_cutoff: float = 3.0, d_h_a_angle_cutoff: float = 150.0, interval: int = 1,) -> None:
def write_salt_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                      sel: str = 'all', hydrogens_sel: str = None, acceptors_sel: str = None,
                      d_a_cutoff: float = 3.0, d_h_a_angle_cutoff: float = 150.0, interval: int = 1,) -> None:
    """Writes hydrogen bonds to `out` file"""
    
    u = mda.Universe(psf, traj)
    hbonds = HBA(universe=u)

    if hydrogens_sel is None or hydrogens_sel.lower() == 'none':
        hbonds.hydrogens_sel = hbonds.guess_hydrogens(sel, max_mass=3.3)
    else:
        hbonds.hydrogens_sel = "({}) and ({})".format(sel, hydrogens_sel)
    if acceptors_sel is None or acceptors_sel.lower() == 'none':
        hbonds.acceptors_sel = hbonds.guess_acceptors(sel)
    else:
        hbonds.acceptors_sel = "({}) and ({})".format(sel, acceptors_sel)
    #if donors_sel is None or donors_sel.lower() == 'none':
    #    hbonds.donors_sel = hbonds.guess_donors(sel)
    #else:
    #    hbonds.donors_sel = "({}) and ({})".format(sel, donors_sel)

    hbonds.d_a_cutoff = d_a_cutoff
    hbonds.d_h_a_angle_cutoff = d_h_a_angle_cutoff

    hbonds.run(step=interval)

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        for entry in hbonds.results['hbonds']:
            #output = "{:.0f} {:.0f} {:.0f} {:.0f} {:f} {:f}".format(*entry)
            output = "{:.0f}".format(entry[0])
            for i in range(1,4):
                atom = u.atoms[int(entry[i])]
                output += ' ' + atom.segid + '_' + atom.resname + str(atom.resid) + '_' + atom.name
            output += " {:f} {:f}".format(*entry[4:])
            print(output, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--d_a_cutoff', type=float, metavar='N', default='3.0',
                        help="Distance cutoff for hydrogen bonds. This cutoff refers to the D-A distance.")
    parser.add_argument('--d_h_a_angle_cutoff', type=float, metavar='N', default='150.0',
                        help="D-H-A angle cutoff for hydrogen bonds (degrees).")
    parser.add_argument('--sel', metavar='selection', default='all',
                        help="Restrict the search to only those atoms")
    parser.add_argument('--hydrogens_sel', metavar='selection', default=None,
                        help="Atom selection for hydrogens. If None, then will be identified via charge and mass.")
    parser.add_argument('--acceptors_sel', metavar='selection', default=None,
                        help="Atom selection for acceptors. If None, then will be identified via charge.")
    ## Only use if the universe topology does not contain bonding information,
    ## otherwise donor-hydrogen pairs may be incorrectly assigned.
    #parser.add_argument('--donors_sel', metavar='selection', default=None,
    #                    help="Atom selection for donors. If None, then will be identified via the topology.")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_salt_bridge(**settings)


if __name__ == '__main__':
    main()

