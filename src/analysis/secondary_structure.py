import argparse
import os
import random
import string
import sys

import MDAnalysis as mda  # type: ignore

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import exec_name

ANALYSIS_NAME = 'secondary_structure'
aa1to3 = {'A': 'ALA', 'V': 'VAL', 'I': 'ILE', 'L': 'LEU',
          'M': 'MET', 'F': 'PHE', 'Y': 'TYR', 'W': 'TRP',
          'S': 'SER', 'T': 'THR', 'N': 'ASN', 'Q': 'GLN',
          'R': 'ARG', 'K': 'LYS', 'H': 'HIS',
          'D': 'ASP', 'E': 'GLU',
          'C': 'CYS', 'G': 'GLY', 'P': 'PRO'}


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "#residue secondary_structures_in_each_frame"

    print(header_str, file=outfile)

    return header_str


def write_salt_bridge(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                      sel: str = 'protein', interval: int = 1,
                      mkdssp_path: str | None = None) -> None:
    """Writes secondary structure to `out` file"""

    if mkdssp_path is None:
        # raises a ValueError if not in PATH
        mkdssp_path = exec_name('mkdssp')

    residue_list: list[str] = []
    ss_all_frames = []
    step_num = 0
    random_tmp_filename = 'tmp' + \
        ''.join(random.choices(string.ascii_letters+string.digits, k=7))

    for traj_file in traj:
        u = mda.Universe(psf, traj_file)
        u.add_TopologyAttr('elements')
        use_segid = False
        if not hasattr(u.atoms, 'chainIDs'):
            use_segid = True
            u.add_TopologyAttr('chainIDs')

        all_atoms = u.select_atoms(sel)
        # Add the column of element that is required by DSSP
        all_atoms.atoms.elements = [name[0] for name in all_atoms.atoms.names]
        # Fix the resnames that are not recognized by DSSP
        resnames = all_atoms.residues.resnames
        for i in range(len(resnames)):
            if resnames[i] in ['HSD', 'HSE', 'HSP']:  # Charmm FF
                resnames[i] = 'HIS'
            if resnames[i] in ['HID', 'HIE', 'HIP']:  # Amber FF
                resnames[i] = 'HIS'
        all_atoms.residues.resnames = resnames
        # Fix the names of C-terminal carboxyl oxygens that are not recognized by DSSP
        atom_names = all_atoms.atoms.names
        for i in range(len(atom_names)):
            # Charmm FF
            if i < len(atom_names) - 1 and atom_names[i] == 'OT1' and atom_names[i+1] == 'OT2':
                atom_names[i] = 'O'
                atom_names[i+1] = 'OXT'
            elif i > 0 and atom_names[i] == 'OT1' and atom_names[i-1] == 'OT2':
                atom_names[i] = 'O'
                atom_names[i-1] = 'OXT'
        all_atoms.atoms.names = atom_names
        if use_segid:
            # If the input system does not has chainIDs, add chainID (A to Z then a to z) to each
            # segment. Could have issue if the system has more than 26 segments.
            chainID_to_segid = {}
            for i in range(all_atoms.n_segments):
                chainID = chr(65+i)
                all_atoms.segments[i].atoms.chainIDs = [
                    chainID]*len(all_atoms.segments[i].atoms)
                chainID_to_segid[chainID] = all_atoms.segids[i]
        for ts in u.trajectory:
            if step_num % interval != 0:
                step_num += 1
                continue
            all_atoms.write('{}_{}.pdb'.format(random_tmp_filename, step_num))
            os.system('{} --output-format=dssp {}_{}.pdb {}_{}.dssp'.format(
                DSSP_path, random_tmp_filename, step_num, random_tmp_filename, step_num))

            processing_line = False
            ss_this_frame = []
            residue_list_this_frame: list[str] = []
            for line in open('{}_{}.dssp'.format(random_tmp_filename, step_num)).readlines():
                if line.strip().startswith("#  RESIDUE AA STRUCTURE"):
                    processing_line = True
                    continue
                if not processing_line:
                    continue
                try:
                    resid = int(line[6:10])
                    chainID = line[11]
                    aa = line[13]
                    ss = line[16]
                    if ss == ' ':
                        ss = '-'
                except:
                    continue
                if len(residue_list) == 0:
                    if use_segid:
                        residue_list_this_frame.append('{}_{}_{}'.format(
                            chainID_to_segid[chainID], aa1to3.get(aa, 'XXX'), resid))
                    else:
                        residue_list_this_frame.append(
                            '{}_{}_{}'.format(chainID, aa, resid))
                ss_this_frame.append(ss)
            ss_all_frames.append(ss_this_frame)
            if len(residue_list) == 0:
                residue_list = residue_list_this_frame
            os.remove('{}_{}.pdb'.format(random_tmp_filename, step_num))
            os.remove('{}_{}.dssp'.format(random_tmp_filename, step_num))
            step_num += 1

    with sta.resolve_file(out, 'w') as outfile:
        print("#DSSP Code  Description", file=outfile)
        print("#    H      Alpha-helix", file=outfile)
        print("#    B      Beta-bridge (residue in isolated beta-bridge)", file=outfile)
        print("#    E      Strand (extended strand, participates in beta-ladder)", file=outfile)
        print("#    G      Helix_3 (3_10-helix)", file=outfile)
        print("#    I      Helix_5 (pi-helix)", file=outfile)
        print("#    P      Helix_PPII (Kappa-helix (poly-proline II helix))", file=outfile)
        print("#    T      Turn (hydrogen-bonded turn)", file=outfile)
        print("#    S      Bend", file=outfile)
        print("#    -      Loop", file=outfile)
        print("#", file=outfile)
        header(outfile)
        for i in range(len(residue_list)):
            output = residue_list[i] + ' '
            for j in range(len(ss_all_frames)):
                output += ss_all_frames[j][i]
            print(output, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    mkdssp_default = exec_name('mkdssp', error=False)

    def mkdssp_or_raise(value: str):
        if value:
            return exec_name(value)

        if mkdssp_default:
            return mkdssp_default

        raise argparse.ArgumentTypeError("No mkdssp in PATH")

    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='protein', default='all',
                        help="Restrict the calculation to only those atoms")

    if mkdssp_default:
        parser.add_argument('--mkdssp-path', metavar='PATH',
                            type=mkdssp_or_raise, default=mkdssp_default,
                            help="Path to mkdssp executable. Default: try to find in PATH.")
    else:
        print("Warning: No mkdssp in PATH", file=sys.stderr)
        parser.add_argument('--mkdssp-path', metavar='PATH',
                            type=mkdssp_or_raise, required=True,
                            help="Path to mkdssp executable. "
                                 "REQUIRED because mkdssp is not in PATH.")

    return parser


def main(settings: dict | None = None) -> None:
    global DSSP_path

    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_salt_bridge(**settings)


if __name__ == '__main__':
    main()
