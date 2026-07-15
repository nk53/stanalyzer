import argparse
import typing as t

import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float
from stanalyzer.workers.contacts_worker import contacts_worker
from stanalyzer.runtime import (
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    chunk_frames,
)

ANALYSIS_NAME = 'contacts'

CPair: t.TypeAlias = tuple[str, int, str, int]
CDict: t.TypeAlias = dict[CPair, int]

def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    header_str = "Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Frequency"

    if not np_formatted:
        header_str = "# " + header_str

    if outfile is not None:
        print(header_str, file=outfile)

    return header_str

def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                   out: sta.FileRef, contact_threshold: float = 5.0,
                   interval: int = 1, debug: bool = False,
                   workers: int | None = None) -> None:
    """
    Calculate residue-residue contact frequencies.

    A contact is recorded when the distance between the centers of mass
    of two residues is below `contact_threshold`.

    Optimized version:
    - Computes residue COMs once per frame.
    - Uses MDAnalysis self_capped_distance() to find only nearby pairs.
    - Avoids O(N^2) Python distance checks for every residue pair.
    """

    if contact_threshold <= 0:
        raise ValueError(
            "contact_threshold must be a positive number, "
            f"not '{contact_threshold}'"
        )

    if interval < 1:
        raise ValueError("interval must be at least 1")

    if workers is not None and workers < 1:
        raise ValueError("workers must be at least 1")

    universe = mda.Universe(psf, traj)

    atoms = universe.select_atoms(sel)
    residues = atoms.residues

    if len(residues) == 0:
        raise ValueError(f"No residues found for selection: {sel}")

    n_frames = len(universe.trajectory)
    context = RuntimeContext.detect_desktop()
    plan = RuntimeScheduler(context).create_plan(
        task_count=n_frames or None,
        n_workers=workers,
    )
    print(
        "Runtime plan: "
        f"strategy={plan.strategy}, "
        f"backend={plan.backend}, workers={plan.n_workers}"
    )

    chunks = chunk_frames(
        n_frames=n_frames,
        n_workers=plan.n_workers,
    )

    tasks = [
        (
            psf,
            traj,
            sel,
            contact_threshold,
            interval,
            debug,
            start,
            stop,
        )
        for start, stop in chunks
    ]

    executor = RuntimeExecutor(plan=plan)
    partial_frequencies = executor.run(
        contacts_worker,
        tasks,
    )

    contact_frequency: CDict = {}
    for partial in partial_frequencies:
        for contact, freq in partial.items():
            contact_frequency[contact] = (
                contact_frequency.get(contact, 0) + freq
            )

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)

        for contact, freq in sorted(contact_frequency.items()):
            print(*contact, freq, file=outfile)

    print(f"Contact frequencies written to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for contact calculation")
    parser.add_argument('--contact-threshold', type=p_float, metavar='N', default='5.0',
                        help="Distance cutoff for calculating the contact frequency.")
    parser.add_argument(
        '--workers',
        type=int,
        metavar='N',
        default=None,
        help=(
            "Maximum number of workers. "
            "Defaults to automatic runtime selection."
        ),
    )
    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()
