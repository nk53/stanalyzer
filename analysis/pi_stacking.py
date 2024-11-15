import argparse
import sys
from typing import Optional
import numpy as np
from operator import itemgetter

import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import normal, angle

ANALYSIS_NAME = 'pi_stacking'


def header(outfile: Optional[sta.FileLike] = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#residue1 residue2 frames"

    print(header_str, file=outfile)

    return header_str


def write_pi_stacking(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                      sel: str = 'all', pi_pi_dist_cutoff: float = 6.0, pi_cation_dist_cutoff: float = 6.0,
                      interval: int = 1,) -> None:
    """Writes pi stacking to `out` file"""
    
    residue2atoms = {'PHE': 'name CG CD* CE* CZ',
                     'TYR': 'name CG CD* CE* CZ',
                     'TRP': 'name CD2 CE2 CZ2 CH2 CZ3 CE3',
                     'HIS': 'name CG ND1 CE1 NE2 CD2'}
    #CHARMM FF
    for residue in ['HSD', 'HSE', 'HSP']:
        residue2atoms[residue] = residue2atoms['HIS']
    #AMBER FF
    for residue in ['HID', 'HIE', 'HIP']:
        residue2atoms[residue] = residue2atoms['HIS']
        residue2atoms['C'+residue] = residue2atoms['HIS']
        residue2atoms['N'+residue] = residue2atoms['HIS']
    
    # the angle between ring normal and vector(ring center to cation) is within 0-30 degrees
    pi_cation_radian_limit = np.pi / 6
    pi_cation_radian_limit2 = np.pi - pi_cation_radian_limit
 
    pi_stacking = {}
    step_num = 0
    
    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
        all_atoms = u.select_atoms(sel)
        
        pi_rings = []
        for residue in all_atoms.residues:
            if residue.resname in residue2atoms:
                pi_ring_atoms = residue.atoms.select_atoms(residue2atoms[residue.resname])
                pi_rings.append(pi_ring_atoms)
        if len(pi_rings) == 0:
            print('Unable to find aromatic residues!')
            return
        cations = all_atoms.select_atoms("(resname ARG LYS CARG CLYS NARG NLYS) and (name NE NH* NZ)")
        if len(pi_rings) + len(cations) < 2:
            print('The total number of aromatic residues and positively charge residues is less than 2!')
            return

        for ts in u.trajectory:
            if step_num % interval != 0:
                step_num += 1
                continue
            pi_ring_centers = []
            pi_ring_normals = []
            for pi_ring in pi_rings:
                pi_ring_center = pi_ring.positions.mean(axis=0)
                pi_ring_centers.append(pi_ring_center)
                v1 = pi_ring.positions[0] - pi_ring_center
                v2 = pi_ring.positions[1] - pi_ring_center
                pi_ring_normals.append(normal(v1, v2))
            pi_ring_centers = np.array(pi_ring_centers)

            pairs, distances = capped_distance(pi_ring_centers, pi_ring_centers,
                                               max_cutoff=pi_pi_dist_cutoff, min_cutoff=1.0, return_distances=True)
            for idx1, idx2 in pairs:
                atom1 = pi_rings[idx1].atoms[0]
                residue1 = atom1.segid + '_' + atom1.resname + '_' + str(atom1.resid)
                atom2 = pi_rings[idx2].atoms[0]
                residue2 = atom2.segid + '_' + atom2.resname + '_' + str(atom2.resid)
                key = tuple(sorted([residue1, residue2]))
                if key not in pi_stacking:
                    pi_stacking[key] = set([step_num])
                else:
                    pi_stacking[key].add(step_num)
            
            pairs2, distances2 = capped_distance(pi_ring_centers, cations,
                                                 max_cutoff=pi_cation_dist_cutoff, min_cutoff=1.0, return_distances=True)
            for idx1, idx2 in pairs2:
                atom1 = pi_rings[idx1].atoms[0]
                residue1 = atom1.segid + '_' + atom1.resname + '_' + str(atom1.resid)
                atom2 = cations[idx2]
                residue2 = atom2.segid + '_' + atom2.resname + '_' + str(atom2.resid)
                # check if the angle between ring normal and vector(ring center to cation) is within 0-30 degrees
                v_center2cation = cations.positions[idx2] - pi_ring_centers[idx1]
                radian = angle(pi_ring_normals[idx1], v_center2cation)
                if (radian >= -pi_cation_radian_limit and radian <= pi_cation_radian_limit
                    or radian >= pi_cation_radian_limit2 or radian <= -pi_cation_radian_limit2):
                    key = (residue1, residue2)
                    if key not in pi_stacking:
                        pi_stacking[key] = set([step_num])
                    else:
                        pi_stacking[key].add(step_num)
 
            step_num += 1

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        count = []
        for key, value in pi_stacking.items():
            count.append([key, len(value)])
        count = sorted(count, key=itemgetter(1), reverse=True)
        for key, _ in count:
            output = "{} {} ".format(*key) + " ".join([str(frame) for frame in sorted(list(pi_stacking[key]))])
            print(output, file=outfile)



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--pi-pi-dist-cutoff', type=float, metavar='N', default='6.0',
                        help="Distance cutoff between aromatic ring centers")
    parser.add_argument('--pi-cation-dist-cutoff', type=float, metavar='N', default='6.0',
                        help="Distance cutoff between aromatic ring centers and positively charged groups")
    parser.add_argument('--sel', metavar='selection', default='all',
                        help="Restrict the search to only those atoms")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_pi_stacking(**settings)


if __name__ == '__main__':
    main()

