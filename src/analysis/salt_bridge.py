import argparse
import numpy as np
from operator import itemgetter

import MDAnalysis as mda  # type: ignore
from MDAnalysis.analysis import contacts  # type: ignore

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

ANALYSIS_NAME = 'salt_bridge'


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#residue1 residue2 frames"

    print(header_str, file=outfile)

    return header_str


def write_salt_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                      positive_sel: str = 'all', negative_sel: str = 'all',
                      positive_def: str | None = None, negative_def: str | None = None,
                      dist_cutoff: float = 4.5, interval: int = 1,) -> None:
    """Writes salt bridge to `out` file"""

    salt_bridges = {}
    step_num = 0

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
        if negative_def is None or negative_def.lower() == 'none':
            acidic = u.select_atoms(negative_sel).select_atoms(
                "(resname ASP GLU CASP CGLU NASP NGLU) and (name OE* OD*)")
        else:
            acidic = u.select_atoms(negative_sel).select_atoms(negative_def)
        if positive_def is None or positive_def.lower() == 'none':
            basic = u.select_atoms(positive_sel).select_atoms(
                "(resname ARG LYS CARG CLYS NARG NLYS) and (name NE NH* NZ)")
        else:
            basic = u.select_atoms(positive_sel).select_atoms(positive_def)
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
            indices1, indices2 = np.where(
                contacts.contact_matrix(dist_mat, dist_cutoff))
            for i in range(len(indices1)):
                atom1 = acidic[indices1[i]]
                residue1 = atom1.segid + '_' + \
                    atom1.resname + '_' + str(atom1.resid)
                atom2 = basic[indices2[i]]
                residue2 = atom2.segid + '_' + \
                    atom2.resname + '_' + str(atom2.resid)
                key = (residue1, residue2)
                if key not in salt_bridges:
                    salt_bridges[key] = set([step_num])
                else:
                    salt_bridges[key].add(step_num)
            step_num += 1

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        count: list = []
        for key, value in salt_bridges.items():
            count.append([key, len(value)])
        count = sorted(count, key=itemgetter(1), reverse=True)
        for key, _ in count:
            output = "{} {} ".format(
                *key) + " ".join([str(frame) for frame in sorted(list(salt_bridges[key]))])
            print(output, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--dist-cutoff', type=p_float, metavar='N', default='4.5',
                        help="Distance cutoff between oppositely charged atoms")
    parser.add_argument('--positive-sel', metavar='selection', default='all',
                        help="Restrict the search of positively charge atoms to only those atoms")
    parser.add_argument('--negative-sel', metavar='selection', default='all',
                        help="Restrict the search of negatively charge atoms to only those atoms")
    parser.add_argument('--positive-def', metavar='selection', default=None,
                        help="Definition of positively charged atoms")
    parser.add_argument('--negative-def', metavar='selection', default=None,
                        help="Definition of negatively charged atoms")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_salt_bridge(**settings)


if __name__ == '__main__':
    main()
