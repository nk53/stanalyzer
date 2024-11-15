import argparse
import sys
from typing import Optional

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

ANALYSIS_NAME = 'hbond'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#frame donor hydrogen acceptor distance angle"

    print(header_str, file=outfile)

    return header_str


#def write_hbond(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
#                sel: str = 'all', hydrogens_sel: str = None, acceptors_sel: str = None, donors_sel: str = None,
#                d_a_cutoff: float = 3.0, d_h_a_angle_cutoff: float = 150.0, interval: int = 1,) -> None:
def write_hbond(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                sel: str = 'all', hydrogens_sel: str = None, acceptors_sel: str = None,
                d_a_cutoff: float = 3.0, d_h_a_angle_cutoff: float = 150.0, interval: int = 1,) -> None:
    """Writes hydrogen bonds to `out` file"""
    
    u = mda.Universe(psf, traj)
    hbonds = HBA(universe=u)

    if hydrogens_sel is None or hydrogens_sel.lower() == 'none':
        try:
            hbonds.hydrogens_sel = hbonds.guess_hydrogens(sel, max_mass=3.3)
        except:
            hbonds.hydrogens_sel = "name HN H H1 H2 H3 HW1 HW2 HT1 HT2 HT3 HN1 HN2 \
                                    or (resname ARG CARG NARG and name HE HH11 HH12 HH21 HH22) \
                                    or (resname ASN CASN NASN and name HD21 HD22) \
                                    or (resname CYS CCYS NCYS and name HG HG1) \
                                    or (resname GLN CGLN NGLN and name HE21 HE22) \
                                    or (resname HSD HID CHID NHID and name HD1) \
                                    or (resname HSE HIE CHIE NHIE and name HE2) \
                                    or (resname HSP HIP CHIP NHIP and name HD1, HE2) \
                                    or (resname LYS CLYS NLYS and name HZ1 HZ2 HZ3) \
                                    or (resname SER CSER NSER and name HG HG1) \
                                    or (resname THR CTHR NTHR and name HG1) \
                                    or (resname TRP CTRP NTRP and name HE1) \
                                    or (resname TYR CTYR NTYR and name HH)"
    else:
        hbonds.hydrogens_sel = "({}) and ({})".format(sel, hydrogens_sel)
    if acceptors_sel is None or acceptors_sel.lower() == 'none':
        try:
            hbonds.acceptors_sel = hbonds.guess_acceptors(sel)
        except:
            hbonds.acceptors_sel = "name O OC1 OC2 OH2 OW* OT1 OT2 \
                                    or (resname ASN CASN NASN and name OD1) \
                                    or (resname ASP CASP NASP and name OD1 OD2) \
                                    or (resname GLN CGLN NGLN and name OE1) \
                                    or (resname GLU CGLU NGLU and name OE1 OE2) \
                                    or (resname HSD HID CHID NHID and name NE2) \
                                    or (resname HSE HIE CHIE NHIE and name ND1) \
                                    or (resname SER CSER NSER and name OG) \
                                    or (resname THR CTHR NTHR and name OG1) \
                                    or (resname TYR CTYR NTYR and name OH)"
            #hbonds.acceptors_sel += "or (resname MET CMET NMET and name SD)"
    else:
        hbonds.acceptors_sel = "({}) and ({})".format(sel, acceptors_sel)
    #if donors_sel is None or donors_sel.lower() == 'none':
    #    try:
    #       hbonds.donors_sel = hbonds.guess_donors(sel)
    #    except:
    #       hbonds.donors_sel = "name N OH2 OW* \
    #                            or (resname ARG CARG NARG and name NE NH1 NH2) \
    #                            or (resname ASN CASN NASN and name ND2) \
    #                            or (resname CYS CCYS NCYS and name SG) \
    #                            or (resname GLN CGLN NGLN and name NE2) \
    #                            or (resname HSD HID CHID NHID and name ND1) \
    #                            or (resname HSE HIE CHIE NHIE and name NE2) \
    #                            or (resname HSP HIP CHIP NHIP and name ND1, NE2) \
    #                            or (resname LYS CLYS NLYS and name NZ) \
    #                            or (resname SER CSER NSER and name OG) \
    #                            or (resname THR CTHR NTHR and name OG1) \
    #                            or (resname TRP CTRP NTRP and name NE1) \
    #                            or (resname TYR CTYR NTYR and name OH)"
    #else:
    #    hbonds.donors_sel = "({}) and ({})".format(sel, donors_sel)

    hbonds.d_a_cutoff = d_a_cutoff
    hbonds.d_h_a_angle_cutoff = d_h_a_angle_cutoff
    #hbonds.d_h_cutoff=1.2

    hbonds.run(step=interval)

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        for entry in hbonds.results['hbonds']:
            #skip if both donor and acceptor are solvent and/or ions
            if ((u.atoms[int(entry[1])].segid == 'SOLV' or u.atoms[int(entry[1])].segid == 'IONS')
                and (u.atoms[int(entry[3])].segid == 'SOLV' or u.atoms[int(entry[3])].segid == 'IONS')):
                continue
            output = "{:.0f}".format(entry[0])
            for i in range(1,4):
                atom = u.atoms[int(entry[i])]
                output += ' ' + atom.segid + '_' + atom.resname + '_' + str(atom.resid) + '_' + atom.name
            output += " {:f} {:f}".format(*entry[4:])
            print(output, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--d-a-cutoff', type=float, metavar='N', default='3.0',
                        help="Distance cutoff for hydrogen bonds. This cutoff refers to the D-A distance.")
    parser.add_argument('--d-h-a-angle-cutoff', type=float, metavar='N', default='150.0',
                        help="D-H-A angle cutoff for hydrogen bonds (degrees).")
    parser.add_argument('--sel', metavar='selection', default='all',
                        help="Restrict the search to only those atoms")
    parser.add_argument('--hydrogens-sel', metavar='selection', default=None,
                        help="Atom selection for hydrogens. If None, then will be identified via charge and mass.")
    parser.add_argument('--acceptors-sel', metavar='selection', default=None,
                        help="Atom selection for acceptors. If None, then will be identified via charge.")
    ## Only use if the universe topology does not contain bonding information,
    ## otherwise donor-hydrogen pairs may be incorrectly assigned.
    #parser.add_argument('--donors-sel', metavar='selection', default=None,
    #                    help="Atom selection for donors. If None, then will be identified via the topology.")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_hbond(**settings)


if __name__ == '__main__':
    main()

