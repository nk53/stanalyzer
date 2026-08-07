import argparse
import re
import typing as t
from collections import defaultdict

import MDAnalysis as mda
import numpy as np

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.stanalyzer import writable_outfile


ANALYSIS_NAME = "bond_statistics"


# ==========================================================
# Types
# ==========================================================

class Bond(t.NamedTuple):
    group1: str
    group2: str


class Angle(t.NamedTuple):
    group1: str
    group2: str
    group3: str


class Dihedral(t.NamedTuple):
    group1: str
    group2: str
    group3: str
    group4: str


OptFileLike: t.TypeAlias = sta.FileRef | None
IndexDict: t.TypeAlias = dict[str, list[int]]
CentroidType: t.TypeAlias = t.Literal["com", "cog"]

BondLengthStats: t.TypeAlias = dict[Bond, list[float]]
BondAngleStats: t.TypeAlias = dict[Angle, list[float]]
BondDihedralStats: t.TypeAlias = dict[Dihedral, list[float]]


class BondParams(t.TypedDict, total=False):
    index: IndexDict
    bonds: list[Bond]
    angles: list[Angle]
    dihedrals: list[Dihedral]


class BondStats(t.TypedDict, total=False):
    bond_lengths: BondLengthStats
    bond_angles: BondAngleStats
    bond_dihedrals: BondDihedralStats


# ==========================================================
# Geometry helpers
# ==========================================================

def calculate_centroid(
    atom_group: mda.AtomGroup,
    method: CentroidType,
) -> np.ndarray:
    """
    Calculate the centroid using MDAnalysis methods so numerical
    behavior remains identical to bond_statistics_v1.
    """

    if method == "cog":
        return atom_group.center_of_geometry()

    if method == "com":
        return atom_group.center_of_mass()

    raise ValueError(
        f"Unknown centroid method '{method}'. "
        "Expected 'cog' or 'com'."
    )


def calculate_bond_length(
    pos1: np.ndarray,
    pos2: np.ndarray,
) -> float:
    """Calculate the distance between two group centroids."""

    return float(np.linalg.norm(pos1 - pos2))


def calculate_bond_angle(
    pos1: np.ndarray,
    pos2: np.ndarray,
    pos3: np.ndarray,
) -> float:
    """
    Calculate the angle formed by three group centroids.

    pos2 is the central point.
    """

    ba = pos1 - pos2
    bc = pos3 - pos2

    norm_ba = np.linalg.norm(ba)
    norm_bc = np.linalg.norm(bc)
    denominator = norm_ba * norm_bc

    if denominator == 0:
        return float("nan")

    cosine_angle = np.dot(ba, bc) / denominator

    angle_rad = np.arccos(
        np.clip(
            cosine_angle,
            -1.0,
            1.0,
        )
    )

    return float(np.degrees(angle_rad))

def calculate_bond_dihedral(
    pos1: np.ndarray,
    pos2: np.ndarray,
    pos3: np.ndarray,
    pos4: np.ndarray,
) -> float:
    """
    Calculate the signed dihedral using the original STAnalyzer
    sign convention.
    """

    ab = pos2 - pos1
    cb = pos3 - pos2
    dc = pos4 - pos3

    normal1 = np.cross(ab, cb)
    normal2 = np.cross(cb, dc)

    norm1 = np.linalg.norm(normal1)
    norm2 = np.linalg.norm(normal2)

    if norm1 == 0 or norm2 == 0:
        return float("nan")

    normal1 = normal1 / norm1
    normal2 = normal2 / norm2

    cosine_phi = np.dot(normal1, normal2)

    phi_rad = np.arccos(
        np.clip(
            cosine_phi,
            -1.0,
            1.0,
        )
    )

    sign = np.sign(
        np.dot(
            np.cross(normal1, normal2),
            cb,
        )
    )

    phi_rad *= sign

    return float(np.degrees(phi_rad))


# ==========================================================
# Input parsing
# ==========================================================

def read_index_file(
    index_file: sta.FileRef,
) -> BondParams:
    """
    Read atom groups, bonds, angles, and dihedrals from an index file.
    """

    def groups_or_error(
        n_expected: int,
        line: str,
    ) -> list[str]:
        groups = line.split()

        if len(groups) != n_expected:
            raise ValueError(
                f"Invalid {section.upper()} entry. "
                f"Expected {n_expected} items, but got {len(groups)}"
            )

        return groups

    groups_indices: IndexDict = defaultdict(list)
    bonds: list[Bond] = []
    angles: list[Angle] = []
    dihedrals: list[Dihedral] = []

    sections = {
        "[INDEX]",
        "[BONDS]",
        "[ANGLES]",
        "[DIHEDRALS]",
    }

    section = ""
    current_group = ""

    with sta.resolve_file(index_file) as infile:
        for raw_line in infile:
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            if line in sections:
                section = line.strip("[]").lower()
                current_group = ""
                continue

            match section:
                case "index":
                    if line.startswith("[") and line.endswith("]"):
                        current_group = line.strip("[]").strip()
                    elif current_group:
                        groups_indices[current_group].extend(
                            map(
                                int,
                                line.split(),
                            )
                        )
                    else:
                        raise ValueError(
                            "Atom indices were found before an INDEX group name"
                        )

                case "bonds":
                    groups = groups_or_error(
                        2,
                        line,
                    )
                    bonds.append(
                        Bond(*groups)
                    )

                case "angles":
                    groups = groups_or_error(
                        3,
                        line,
                    )
                    angles.append(
                        Angle(*groups)
                    )

                case "dihedrals":
                    groups = groups_or_error(
                        4,
                        line,
                    )
                    dihedrals.append(
                        Dihedral(*groups)
                    )

                case _:
                    raise ValueError(
                        f"Content found outside a recognized section: {line}"
                    )

    if not groups_indices:
        raise ValueError("No atom groups were found in the index file")

    parameters: BondParams = {
        "index": dict(groups_indices),
    }

    if bonds:
        parameters["bonds"] = bonds

    if angles:
        parameters["angles"] = angles

    if dihedrals:
        parameters["dihedrals"] = dihedrals

    validate_parameter_groups(parameters)

    return parameters


def convert_atom_groups(
    input_str: str,
) -> BondParams:
    """
    Convert compact CLI syntax into group and geometry definitions.

    Examples
    --------
    Bond:
        (1,2,3)(4,5,6)

    Angle:
        (1,2)(3,4)(5,6)

    Dihedral:
        (1)(2)(3)(4)
    """

    raw_groups = re.findall(
        r"\((.*?)\)",
        input_str,
    )

    if not raw_groups:
        raise ValueError(
            "No atom groups were found. "
            "Expected syntax such as '(1,2,3)(4,5,6)'."
        )

    groups_indices: IndexDict = {}

    for index, raw_group in enumerate(
        raw_groups,
        start=1,
    ):
        group_name = f"G{index}"

        try:
            atom_indices = [
                int(value.strip())
                for value in raw_group.split(",")
                if value.strip()
            ]
        except ValueError as error:
            raise ValueError(
                f"Invalid atom index in group {group_name}: {raw_group}"
            ) from error

        if not atom_indices:
            raise ValueError(
                f"Group {group_name} contains no atom indices"
            )

        groups_indices[group_name] = atom_indices

    parameters: BondParams = {
        "index": groups_indices,
    }

    match len(raw_groups):
        case 2:
            parameters["bonds"] = [
                Bond(
                    "G1",
                    "G2",
                )
            ]

        case 3:
            parameters["angles"] = [
                Angle(
                    "G1",
                    "G2",
                    "G3",
                )
            ]

        case 4:
            parameters["dihedrals"] = [
                Dihedral(
                    "G1",
                    "G2",
                    "G3",
                    "G4",
                )
            ]

        case n_groups:
            raise ValueError(
                "Invalid number of groups: "
                f"{n_groups}. Expected 2, 3, or 4."
            )

    return parameters


def validate_parameter_groups(
    parameters: BondParams,
) -> None:
    """
    Verify that all bond, angle, and dihedral definitions reference
    existing atom groups.
    """

    group_names = set(
        parameters.get(
            "index",
            {},
        )
    )

    if not group_names:
        raise ValueError("No atom groups were defined")

    definitions: list[tuple[str, tuple[str, ...]]] = []

    for bond in parameters.get("bonds", []):
        definitions.append(
            (
                "bond",
                tuple(bond),
            )
        )

    for angle in parameters.get("angles", []):
        definitions.append(
            (
                "angle",
                tuple(angle),
            )
        )

    for dihedral in parameters.get("dihedrals", []):
        definitions.append(
            (
                "dihedral",
                tuple(dihedral),
            )
        )

    for definition_type, referenced_groups in definitions:
        missing = [
            group
            for group in referenced_groups
            if group not in group_names
        ]

        if missing:
            raise ValueError(
                f"{definition_type.capitalize()} references undefined "
                f"group(s): {', '.join(missing)}"
            )


def process_bond_parameters(
    filename_or_str: sta.FileRef,
    index: bool = False,
) -> BondParams:
    """
    Parse either an index file or compact atom-group string.
    """

    if index:
        parameters = read_index_file(
            filename_or_str
        )
    elif isinstance(filename_or_str, str) and filename_or_str:
        parameters = convert_atom_groups(
            filename_or_str
        )
    else:
        raise ValueError(
            "Need either an index file or an atom-group string"
        )

    validate_parameter_groups(
        parameters
    )

    return parameters


# ==========================================================
# Streaming analysis
# ==========================================================

def build_atom_groups(
    universe: mda.Universe,
    groups: IndexDict,
) -> dict[str, mda.AtomGroup]:
    """
    Construct each AtomGroup once before trajectory iteration.

    Input atom indices are one-based, so they are converted to zero-based
    indices exactly once here.
    """

    atom_groups: dict[str, mda.AtomGroup] = {}
    n_atoms = len(universe.atoms)

    for group_name, indices in groups.items():
        if not indices:
            raise ValueError(
                f"Atom group '{group_name}' contains no indices"
            )

        zero_based = np.asarray(
            indices,
            dtype=np.int64,
        ) - 1

        if np.any(zero_based < 0) or np.any(zero_based >= n_atoms):
            raise IndexError(
                f"Atom group '{group_name}' contains indices outside "
                f"the valid one-based range 1-{n_atoms}"
            )

        atom_groups[group_name] = universe.atoms[
            zero_based
        ]

    return atom_groups


def analyze_bond_parameters_from_process(
    universe: mda.Universe,
    index_file: OptFileLike = None,
    atom_groups: str = "",
    centroid: CentroidType = "cog",
) -> BondStats:
    """
    Analyze bonds, angles, and dihedrals in a single trajectory pass.

    Phase 1 design:
    - Parse geometry definitions once.
    - Construct MDAnalysis AtomGroups once.
    - Cache masses once for COM calculations.
    - Calculate group centroids once per frame.
    - Calculate requested metrics immediately.
    - Do not store the complete centroid trajectory.
    """

    filename_or_str = (
        index_file
        if index_file is not None
        else atom_groups
    )

    parameters = process_bond_parameters(
        filename_or_str,
        index=index_file is not None,
    )

    groups = parameters.get(
        "index",
        {},
    )
    bonds = parameters.get(
        "bonds",
        [],
    )
    angles = parameters.get(
        "angles",
        [],
    )
    dihedrals = parameters.get(
        "dihedrals",
        [],
    )

    group_atomgroups = build_atom_groups(
        universe,
        groups,
    )

    bond_lengths: BondLengthStats = defaultdict(list)
    bond_angles: BondAngleStats = defaultdict(list)
    bond_dihedrals: BondDihedralStats = defaultdict(list)

    for _ in universe.trajectory:
        # Calculate every group centroid once for this frame.
        frame_centroids = {
            name: calculate_centroid(
                atom_group,
                method=centroid,
            )
            for name, atom_group in group_atomgroups.items()
        }

        for bond in bonds:
            value = calculate_bond_length(
                frame_centroids[bond.group1],
                frame_centroids[bond.group2],
            )
            bond_lengths[bond].append(
                value
            )

        for angle_definition in angles:
            value = calculate_bond_angle(
                frame_centroids[angle_definition.group1],
                frame_centroids[angle_definition.group2],
                frame_centroids[angle_definition.group3],
            )
            bond_angles[angle_definition].append(
                value
            )

        for dihedral in dihedrals:
            value = calculate_bond_dihedral(
                frame_centroids[dihedral.group1],
                frame_centroids[dihedral.group2],
                frame_centroids[dihedral.group3],
                frame_centroids[dihedral.group4],
            )
            bond_dihedrals[dihedral].append(
                value
            )

    results: BondStats = {}

    if bond_lengths:
        results["bond_lengths"] = dict(
            bond_lengths
        )

    if bond_angles:
        results["bond_angles"] = dict(
            bond_angles
        )

    if bond_dihedrals:
        results["bond_dihedrals"] = dict(
            bond_dihedrals
        )

    return results


# ==========================================================
# Output
# ==========================================================

def write_bond_lengths_to_dat(
    outfile: sta.FileRef,
    bond_lengths: BondLengthStats,
) -> None:
    """Write bond lengths to a data file."""

    if not bond_lengths:
        return

    with sta.resolve_file(
        outfile,
        "w",
    ) as output:
        for bond, lengths in bond_lengths.items():
            label = (
                f"[{bond.group1}_{bond.group2}]"
            )

            output.write(
                f"@Bond Length (Angstrom){label}\n"
            )

            for length in lengths:
                output.write(
                    f"{length:.4f}\n"
                )


def write_bond_angles_to_dat(
    outfile: sta.FileRef,
    bond_angles: BondAngleStats,
) -> None:
    """Write bond angles to a data file."""

    if not bond_angles:
        return

    with sta.resolve_file(
        outfile,
        "w",
    ) as output:
        for angle_definition, values in bond_angles.items():
            label = (
                f"[{angle_definition.group1}_"
                f"{angle_definition.group2}_"
                f"{angle_definition.group3}]"
            )

            output.write(
                f"@Bond Angle (Degrees){label}\n"
            )

            for value in values:
                output.write(
                    f"{value:.4f}\n"
                )


def write_bond_dihedrals_to_dat(
    outfile: sta.FileRef,
    bond_dihedrals: BondDihedralStats,
) -> None:
    """Write bond dihedrals to a data file."""

    if not bond_dihedrals:
        return

    with sta.resolve_file(
        outfile,
        "w",
    ) as output:
        for dihedral, values in bond_dihedrals.items():
            label = (
                f"[{dihedral.group1}_"
                f"{dihedral.group2}_"
                f"{dihedral.group3}_"
                f"{dihedral.group4}]"
            )

            output.write(
                f"@Bond Dihedral (Degrees){label}\n"
            )

            for value in values:
                output.write(
                    f"{value:.4f}\n"
                )


def write_files(
    results: BondStats,
    bond_out: sta.FileRef = "bond_lengths.dat",
    angle_out: sta.FileRef = "bond_angles.dat",
    dihedral_out: sta.FileRef = "bond_dihedrals.dat",
) -> None:
    """Write each available analysis result."""

    bond_lengths = results.get(
        "bond_lengths"
    )

    if bond_lengths:
        write_bond_lengths_to_dat(
            bond_out,
            bond_lengths,
        )

    bond_angles = results.get(
        "bond_angles"
    )

    if bond_angles:
        write_bond_angles_to_dat(
            angle_out,
            bond_angles,
        )

    bond_dihedrals = results.get(
        "bond_dihedrals"
    )

    if bond_dihedrals:
        write_bond_dihedrals_to_dat(
            dihedral_out,
            bond_dihedrals,
        )


# ==========================================================
# CLI
# ==========================================================

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=f"stanalyzer {ANALYSIS_NAME}"
    )

    sta.add_project_args(
        parser,
        "psf",
        "traj",
    )

    parser.add_argument(
        "-c",
        "--centroid",
        metavar="OPT",
        default="cog",
        choices=[
            "cog",
            "com",
        ],
        help=(
            "Centroid calculation method. "
            "com: center of mass; cog: center of geometry. "
            "Default: cog."
        ),
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "-a",
        "--atom-groups",
        metavar="GROUPS",
        help=(
            "Atom indices for groups to analyze. "
            "Two groups define a bond, three define an angle, "
            "and four define a dihedral. "
            "Example: (1,2,3)(4,5,6)"
        ),
    )

    group.add_argument(
        "-i",
        "--index-file",
        metavar="FILE",
        type=sta.InputFile,
        help="File containing named atom groups and geometry definitions.",
    )

    parser.add_argument(
        "-bo",
        "--bond-out",
        metavar="FILE",
        type=writable_outfile,
        default="bond_lengths.dat",
        help="Bond-length output file.",
    )

    parser.add_argument(
        "-ao",
        "--angle-out",
        metavar="FILE",
        type=writable_outfile,
        default="bond_angles.dat",
        help="Bond-angle output file.",
    )

    parser.add_argument(
        "-do",
        "--dihedral-out",
        metavar="FILE",
        type=writable_outfile,
        default="bond_dihedrals.dat",
        help="Dihedral-angle output file.",
    )

    return parser


def main(
    settings: dict | None = None,
) -> None:
    if settings is None:
        settings = dict(
            sta.get_settings(
                ANALYSIS_NAME
            )
        )

    psf = settings.pop(
        "psf"
    )
    traj = settings.pop(
        "traj"
    )

    outfiles = {
        "bond_out": settings.pop(
            "bond_out"
        ),
        "angle_out": settings.pop(
            "angle_out"
        ),
        "dihedral_out": settings.pop(
            "dihedral_out"
        ),
    }

    universe = mda.Universe(
        psf,
        traj,
    )

    results = analyze_bond_parameters_from_process(
        universe=universe,
        **settings,
    )

    write_files(
        results,
        **outfiles,
    )


if __name__ == "__main__":
    main()