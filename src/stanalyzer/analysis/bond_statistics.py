import argparse
import re
import typing as t
from collections import defaultdict, namedtuple
from collections.abc import Iterable

import numpy as np
import MDAnalysis as mda
import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.stanalyzer import writable_outfile

if t.TYPE_CHECKING:
    import numpy.typing as npt

ANALYSIS_NAME = 'bond_statistics'

# defining bond, angles and dihedral tuples
Bond = namedtuple('Bond', ['group1', 'group2'])
Angle = namedtuple('Angle', ['group1', 'group2', 'group3'])
Dihedral = namedtuple('Dihedral', ['group1', 'group2', 'group3', 'group4'])

T = t.TypeVar('T')  # for generic type definitions
G = t.TypeVar('G', Bond, Angle, Dihedral)  # G = group

OptFileLike: t.TypeAlias = sta.FileRef | None
Params: t.TypeAlias = list[Bond] | list[Angle] | list[Dihedral]
IndexDict: t.TypeAlias = dict[str, list[int]]
Coors: t.TypeAlias = 'npt.NDArray[np.float64]'
CoorMap: t.TypeAlias = dict[str, list[Coors]]
CentroidType: t.TypeAlias = t.Literal['com', 'cog']
Stats: t.TypeAlias = dict[G, list[np.float64]]


class BondParams(t.TypedDict, total=False):
    index: IndexDict
    bonds: list[Bond]
    angles: list[Angle]
    dihedrals: list[Dihedral]


class BondStats(t.TypedDict, total=False):
    atomgroup_positions: CoorMap
    bond_lengths: Stats[Bond]
    bond_angles: Stats[Angle]
    bond_dihedrals: Stats[Dihedral]


def calculate_atomgroup_positions(
        u: mda.Universe, groups: IndexDict,
        method: CentroidType = "cog") -> CoorMap:
    """Calculates the position of each group of atoms using center of geometry
    or center of mass based on user input."""
    atomgroup_positions: dict[str, list[Coors]] = defaultdict(list)
    centroid: Coors
    for ts in u.trajectory:
        for group, indices in groups.items():
            selected_atoms = u.atoms[np.array(indices) - 1]
            if method == "cog":
                centroid = selected_atoms.center_of_geometry()
            elif method == "com":
                centroid = selected_atoms.center_of_mass()

            atomgroup_positions[group].append(centroid)
    return atomgroup_positions


def calculate_bond_lengths(atomgroup_positions: CoorMap,
                           bonds: Iterable[Bond]) -> Stats[Bond]:
    """Calculates bond lengths for each bond and frame."""
    bond_lengths: Stats[Bond] = defaultdict(list)
    for bond in bonds:
        for pos1, pos2 in zip(atomgroup_positions[bond.group1], atomgroup_positions[bond.group2]):
            length = np.linalg.norm(pos1 - pos2)
            bond_lengths[bond].append(length)
    return bond_lengths


def calculate_bond_angles(atomgroup_positions: CoorMap,
                          angles: Iterable[Angle]) -> Stats[Angle]:
    """Calculates bond angles for each angle and frame."""
    bond_angles: Stats[Angle] = defaultdict(list)
    for angle in angles:
        for pos1, pos2, pos3 in zip(atomgroup_positions[angle.group1],
                                    atomgroup_positions[angle.group2],
                                    atomgroup_positions[angle.group3]):
            ba = pos1 - pos2
            bc = pos3 - pos2
            cosine_angle = np.dot(ba, bc) / \
                (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle_rad = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
            angle_deg = np.degrees(angle_rad)
            bond_angles[angle].append(angle_deg)
    return bond_angles


def calculate_bond_dihedrals(atomgroup_positions: CoorMap,
                             dihedrals: Iterable[Dihedral]) -> Stats[Dihedral]:
    """Calculates dihedral angles for each dihedral and frame."""
    dihedral_angles = defaultdict(list)
    for dihedral in dihedrals:
        for pos1, pos2, pos3, pos4 in zip(atomgroup_positions[dihedral.group1],
                                          atomgroup_positions[dihedral.group2],
                                          atomgroup_positions[dihedral.group3],
                                          atomgroup_positions[dihedral.group4]):
            ab = pos2 - pos1
            cb = pos3 - pos2
            dc = pos4 - pos3

            # Normal vectors to the planes
            normal1 = np.cross(ab, cb)
            normal2 = np.cross(cb, dc)
            # Normalize the normal vectors
            normal1 /= np.linalg.norm(normal1)
            normal2 /= np.linalg.norm(normal2)

            # Cosine of the angle
            cosine_phi = np.dot(normal1, normal2)
            phi_rad = np.arccos(np.clip(cosine_phi, -1.0, 1.0))

            # Calculate the sign of the angle using the direction of cb
            sign = np.sign(np.dot(np.cross(normal1, normal2), cb))
            phi_rad *= sign
            phi_deg = np.degrees(phi_rad)
            dihedral_angles[dihedral].append(phi_deg)
    return dihedral_angles


def read_index_file(index_file: sta.FileRef) -> BondParams:
    """Reads the index file and returns a dictionary of bead groups, bonds,
    angles, and dihedrals."""

    def groups_or_error(n_expected: int, line: str) -> list[str]:
        assert in_section

        groups = line.split()
        n_groups = len(groups)
        if n_groups != n_expected:
            tpl = "Invalid {} entry. Expected {} items, but got {}"
            raise ValueError(tpl.format(section.upper(), n_expected, n_groups))

        return groups

    groups_indices: IndexDict = defaultdict(list)
    bonds: list[Bond] = []
    angles: list[Angle] = []
    dihedrals: list[Dihedral] = []

    sections = '[INDEX]', '[BONDS]', '[ANGLES]', '[DIHEDRALS]'
    in_section = False
    current_group = ''
    section = ''

    with sta.resolve_file(index_file) as f:
        for line in f:
            line = line.strip()

            if line in sections:
                section = line.strip('[]').lower()
                in_section = True
                continue

            match section:
                case 'index':
                    if line.startswith('['):
                        current_group = line.strip('[]').strip()
                    elif line.strip() and current_group:
                        groups_indices[current_group].extend(
                            map(int, line.split()))
                case 'bonds':
                    groups = groups_or_error(2, line)
                    bonds.append(Bond(*groups))
                case 'angles':
                    groups = groups_or_error(3, line)
                    angles.append(Angle(*groups))
                case 'dihedrals':
                    groups = groups_or_error(4, line)
                    dihedrals.append(Dihedral(*groups))

    # Return results, excluding empty lists
    bond_parameters: BondParams = {'index': groups_indices}
    if bonds:
        bond_parameters['bonds'] = bonds
    if angles:
        bond_parameters['angles'] = angles
    if dihedrals:
        bond_parameters['dihedrals'] = dihedrals

    return bond_parameters


def convert_atom_groups(input_str: str) -> BondParams:
    """Extracts information from text input"""
    # Initialize the defaultdict for groups and dictionary for results
    groups_indices: IndexDict = defaultdict(list)

    # Use regex to extract the groups from the string
    groups: list[str] = re.findall(r'\((.*?)\)', input_str)

    # Convert list[str] into groups
    for i, group in enumerate(groups):
        group_name = f'G{i+1}'
        group_list = list(map(int, group.split(',')))
        groups_indices[group_name] = group_list

    bond_parameters: BondParams = {'index': groups_indices}

    # Prepare the bond parameters based on the method
    match len(groups):
        case 2:
            bond_parameters['bonds'] = [Bond(group1='G1', group2='G2')]
        case 3:
            bond_parameters['angles'] = [
                Angle(group1='G1', group2='G2', group3='G3')]
        case 4:
            bond_parameters['dihedrals'] = [
                Dihedral(group1='G1', group2='G2', group3='G3', group4='G4')]
        case n_groups:
            raise ValueError(f"Invalid group count: {n_groups}")

    # Include the index section, as required
    bond_parameters['index'].update(groups_indices)

    # Return results, excluding empty lists (similar to `read_index_file`)
    return bond_parameters


def process_bond_parameters(filename_or_str: sta.FileRef,
                            index: bool = False) -> BondParams:
    """Determines whether to run read_index_file or convert_atom_groups based
    on the index parameter.

    :param filename_or_str: Filename for reading or input string for conversion
    :param method: Method for convert_atom_groups, defaults to 'Bonds' if index is False
    :param index: Boolean to decide which function to run
    :return: Bond parameters from the respective function
    """
    if index:
        return read_index_file(filename_or_str)
    if isinstance(filename_or_str, str) and filename_or_str:
        return convert_atom_groups(filename_or_str)

    raise ValueError("Need either an index file or an atom group str")


def analyze_bond_parameters_from_process(
        universe: mda.Universe, index_file: OptFileLike = None,
        atom_groups: str = '', centroid: CentroidType = 'cog') -> BondStats:

    """Determines whether to run read_index_file or convert_atom_groups based
    on the index parameter, then calculates atom group positions, bond lengths,
    bond angles, and bond dihedrals.

    :param filename_or_str: Filename for reading or input string for conversion
    :param method: Method for convert_atom_groups, defaults to 'Bonds' if index is False
    :param index: Boolean to decide which function to run
    :param universe: MDAnalysis universe object containing trajectory and atom
                     information (required if index is False)
    :param calculation_method: Method to calculate group positions, "cog"
                               (center of geometry) or "com" (center of mass)

    :return: Dictionary with atom group positions, bond lengths, bond angles,
             and dihedrals
    """

    filename_or_str = index_file or atom_groups
    index = index_file is not None

    # Process the bond parameters
    bond_parameters = process_bond_parameters(filename_or_str, index=index)

    # Extract data from the returned dictionary
    groups = bond_parameters.get('index', {})
    bonds = bond_parameters.get('bonds', [])
    angles = bond_parameters.get('angles', [])
    dihedrals = bond_parameters.get('dihedrals', [])

    # Calculate atom group positions (Center of Geometry or Center of Mass)
    atomgroup_positions = calculate_atomgroup_positions(
        universe, groups, method=centroid)

    # Initialize results dictionary
    results: BondStats = {'atomgroup_positions': atomgroup_positions}

    if atomgroup_positions:
        # Calculate bond lengths if bonds are provided
        if bonds:
            bond_lengths = calculate_bond_lengths(atomgroup_positions, bonds)
            results['bond_lengths'] = bond_lengths

        # Calculate bond angles if angles are provided
        if angles:
            bond_angles = calculate_bond_angles(atomgroup_positions, angles)
            results['bond_angles'] = bond_angles

        # Calculate bond dihedrals if dihedrals are provided
        if dihedrals:
            bond_dihedrals = calculate_bond_dihedrals(
                atomgroup_positions, dihedrals)
            results['bond_dihedrals'] = bond_dihedrals

    return results


def write_bond_lengths_to_dat(outfile: sta.FileRef,
                              bond_lengths: Stats[Bond]) -> None:
    """Writes bond lengths to a .dat file if they exist."""
    if bond_lengths:  # Check if there are bond lengths to write
        with sta.resolve_file(outfile, 'w') as f:
            for bond, lengths in bond_lengths.items():
                indices = f"[{''.join(map(str, bond.group1))}_{''.join(map(str, bond.group2))}]"
                f.write(f"@Bond Length (Angstrom){indices}\n")
                for length in lengths:
                    f.write(f"{length:.4f}\n")


def write_bond_angles_to_dat(outfile: sta.FileRef,
                             bond_angles: Stats[Angle]) -> None:
    """Writes bond angles to a .dat file if they exist."""
    if bond_angles:  # Check if there are bond angles to write
        with sta.resolve_file(outfile, 'w') as f:
            for angle, values in bond_angles.items():
                indices = f"[{''.join(map(str, angle.group1))}" \
                          f"_{''.join(map(str, angle.group2))}" \
                          f"_{''.join(map(str, angle.group3))}]"
                f.write(f"@Bond Angle (Degrees){indices}\n")
                for value in values:
                    f.write(f"{value:.4f}\n")


def write_bond_dihedrals_to_dat(outfile: sta.FileRef,
                                bond_dihedrals: Stats[Dihedral]) -> None:
    """Writes bond dihedrals to a .dat file if they exist."""
    if bond_dihedrals:  # Check if there are bond dihedrals to write
        with sta.resolve_file(outfile, 'w') as f:
            for dihedral, values in bond_dihedrals.items():
                indices = f"[{''.join(map(str, dihedral.group1))}" \
                          f"_{''.join(map(str, dihedral.group2))}" \
                          f"_{''.join(map(str, dihedral.group3))}" \
                          f"_{''.join(map(str, dihedral.group4))}]"
                f.write(f"@Bond Dihedral (Degrees){indices}\n")
                for value in values:
                    f.write(f"{value:.4f}\n")


def write_files(results: BondStats,
                bond_out: sta.FileRef = 'bond_lengths.dat',
                angle_out: sta.FileRef = 'bond_angles.dat',
                dihedral_out: sta.FileRef = 'bond_dihedrals.dat') -> None:
    """Write to individual .dat files if data exists"""
    if results.get('bond_lengths'):
        write_bond_lengths_to_dat(
            bond_out, bond_lengths=results['bond_lengths'])

    if results.get('bond_angles'):
        write_bond_angles_to_dat(
            angle_out, bond_angles=results['bond_angles'])

    if results.get('bond_dihedrals'):
        write_bond_dihedrals_to_dat(
            dihedral_out, bond_dihedrals=results['bond_dihedrals'])


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj')
    parser.add_argument('-c', '--centroid', metavar='OPT', default='cog', choices=['cog', 'com'],
                        help="com: center of mass; cog: center of geometry. Default: cog")

    # either -a or -i must be passed, but not both
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--atom-groups', metavar='GROUPS',
                       help="Atom indices for groups to analyze. Number of groups given "
                       "determines analysis type. 2: bond, 3: angle, 4: dihedral. "
                       "Bond example: (1,2,3)(4,5,6).")
    group.add_argument('-i', '--index-file', metavar='FILE', type=sta.InputFile,
                       help="File containing indices to read.")
    parser.add_argument('-bo', '--bond-out', metavar='FILE', type=writable_outfile,
                        default='bond_lengths.dat', help="Location to write bonds.")
    parser.add_argument('-ao', '--angle-out', metavar='FILE', type=writable_outfile,
                        default='bond_angles.dat', help="Location to write bonds.")
    parser.add_argument('-do', '--dihedral-out', metavar='FILE', type=writable_outfile,
                        default='bond_dihedrals.dat', help="Location to write bonds.")

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    u = mda.Universe(settings.pop('psf'), settings.pop('traj'))
    outfiles = {(k := f"{t}_out"): settings.pop(k)
                for t in ('bond', 'angle', 'dihedral')}

    results = analyze_bond_parameters_from_process(universe=u, **settings)
    write_files(results, **outfiles)


if __name__ == "__main__":
    main()
