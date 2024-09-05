import argparse
import sys
from typing import Optional
import numpy as np
from operator import itemgetter

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

ANALYSIS_NAME = 'salt_bridge'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#residue1 residue2 frames"

    print(header_str, file=outfile)

    return header_str


def write_salt_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                      sel: str = 'all', positive_sel: str = None, negative_sel: str = None,
                      dist_cutoff: float = 4.5, interval: int = 1,) -> None:
    """Writes salt bridge to `out` file"""
    
    salt_bridges = {}
    step_num = 0

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
        atoms = u.select_atoms(sel)
        if negative_sel is None:
            acidic = atoms.select_atoms("(resname ASP GLU) and (name OE* OD*)")
        else:
            acidic = atoms.select_atoms(negative_sel)
        if positive_sel is None:
            basic = atoms.select_atoms("(resname ARG LYS) and (name NH* NZ)")
        else:
            basic = atoms.select_atoms(positive_sel)
        if len(acidic) == 0:
            print('Unable to find negatively charged atoms!')
            return
        if len(basic) == 0:
            print('Unable to find positively charged atoms!')
            return
        for ts in u.trajectory:
            if step_num % interval != 0:
                step_num += 1
                continue
            dist_mat = contacts.distance_array(acidic, basic)
            indices1, indices2 = np.where(contacts.contact_matrix(dist_mat, dist_cutoff))
            for i in range(len(indices1)):
                atom1 = acidic[indices1[i]]
                residue1 = atom1.segid + '_' + atom1.resname + str(atom1.resid)
                atom2 = basic[indices2[i]]
                residue2 = atom2.segid + '_' + atom2.resname + str(atom2.resid)
                key = (residue1, residue2)
                if key not in salt_bridges:
                    salt_bridges[key] = set([step_num])
                else:
                    salt_bridges[key].add(step_num)
            step_num += 1

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        count = []
        for key, value in salt_bridges.items():
            count.append([key, len(value)])
        count = sorted(count, key=itemgetter(1), reverse=True)
        for key, _ in count:
            output = "{} {} ".format(*key) + " ".join([str(frame) for frame in sorted(list(salt_bridges[key]))])
            print(output, file=outfile)



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--dist_cutoff', type=float, metavar='N', default='4.5',
                        help="Distance cutoff between oppositely charged atoms")
    parser.add_argument('--sel', metavar='selection', default='all',
                        help="Restrict the search to only those atoms")
    parser.add_argument('--positive_sel', metavar='selection', default=None,
                        help="Atom selection for positively charged atoms")
    parser.add_argument('--negative_sel', metavar='selection', default=None,
                        help="Atom selection for negatively charged atoms")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_salt_bridge(**settings)


if __name__ == '__main__':
    main()

