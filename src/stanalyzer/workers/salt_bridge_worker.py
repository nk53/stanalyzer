from collections import defaultdict

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance


SaltBridge = tuple[str, str]
SaltBridgeFrames = dict[SaltBridge, set[int]]


DEFAULT_NEGATIVE_DEF = (
    "(resname ASP GLU CASP CGLU NASP NGLU) "
    "and (name OE* OD*)"
)

DEFAULT_POSITIVE_DEF = (
    "(resname ARG LYS CARG CLYS NARG NLYS) "
    "and (name NE NH* NZ)"
)


def build_selection(
    outer_selection: str,
    atom_definition: str | None,
    default_definition: str,
) -> str:
    """Combine a user scope with a charged-atom definition."""

    definition = (
        default_definition
        if atom_definition is None
        or atom_definition.lower() == "none"
        else atom_definition
    )

    return f"({outer_selection}) and ({definition})"


def salt_bridge_worker(task) -> SaltBridgeFrames:
    """
    Calculate salt bridges for one frame chunk.

    Parameters
    ----------
    task
        (
            psf,
            traj,
            positive_selection,
            negative_selection,
            dist_cutoff,
            interval,
            start,
            stop,
            debug,
        )

    Returns
    -------
    dict
        Mapping:

            (negative_residue, positive_residue) -> {frame indices}
    """

    (
        psf,
        traj,
        positive_selection,
        negative_selection,
        dist_cutoff,
        interval,
        start,
        stop,
        debug,
    ) = task

    universe = mda.Universe(psf, traj)

    acidic = universe.select_atoms(negative_selection)
    basic = universe.select_atoms(positive_selection)

    if len(acidic) == 0:
        raise ValueError(
            "Unable to find negatively charged atoms with selection: "
            f"{negative_selection}"
        )

    if len(basic) == 0:
        raise ValueError(
            "Unable to find positively charged atoms with selection: "
            f"{positive_selection}"
        )

    # Cache residue labels once. Avoid repeated segid/resname/resid access
    # for every pair in every frame.
    acidic_labels = [
        f"{atom.segid}_{atom.resname}_{atom.resid}"
        for atom in acidic
    ]

    basic_labels = [
        f"{atom.segid}_{atom.resname}_{atom.resid}"
        for atom in basic
    ]

    salt_bridges: dict[SaltBridge, set[int]] = defaultdict(set)

    # Keep interval aligned to global trajectory indices.
    first_frame = start + (-start % interval)

    for frame_idx in range(first_frame, stop, interval):
        ts = universe.trajectory[frame_idx]

        pairs = capped_distance(
            acidic.positions,
            basic.positions,
            max_cutoff=dist_cutoff,
            box=ts.dimensions,
            return_distances=False,
        )

        # Multiple charged atoms may identify the same residue-residue
        # bridge in one frame. Count the residue pair only once per frame.
        frame_bridges: set[SaltBridge] = set()

        for acidic_idx, basic_idx in pairs:
            frame_bridges.add(
                (
                    acidic_labels[int(acidic_idx)],
                    basic_labels[int(basic_idx)],
                )
            )

        for bridge in frame_bridges:
            salt_bridges[bridge].add(frame_idx)

        if debug:
            print(
                f"frame={frame_idx} "
                f"atom_pairs={len(pairs)} "
                f"residue_bridges={len(frame_bridges)}"
            )

    return dict(salt_bridges)