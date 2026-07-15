import argparse
import time
from collections import defaultdict

import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta

from stanalyzer.workers.pi_stacking_worker import (
    CATION_SELECTION,
    RESIDUE_TO_RING_ATOMS,
    pi_stacking_worker,
)
from stanalyzer.runtime import (
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    chunk_frames,
)


ANALYSIS_NAME = "pi_stacking"

PiPair = tuple[str, str]
PiEvents = dict[PiPair, set[int]]


def header(
    outfile: sta.FileLike | None = None,
) -> str:
    header_str = "#residue1 residue2 frames"

    if outfile is not None:
        print(header_str, file=outfile)

    return header_str


def merge_events(
    partial_results: list[PiEvents],
) -> PiEvents:
    merged: dict[PiPair, set[int]] = defaultdict(set)

    for partial in partial_results:
        for pair, frames in partial.items():
            merged[pair].update(frames)

    return dict(merged)


def write_pi_stacking(
    psf: sta.FileRef,
    traj: sta.FileRefList,
    out: sta.FileRef,
    sel: str = "all",
    pi_pi_dist_cutoff: float = 6.0,
    pi_cation_dist_cutoff: float = 6.0,
    interval: int = 1,
    workers: int | None = None,
    debug: bool = False,
) -> None:
    if pi_pi_dist_cutoff <= 0:
        raise ValueError(
            "pi_pi_dist_cutoff must be positive"
        )

    if pi_cation_dist_cutoff <= 0:
        raise ValueError(
            "pi_cation_dist_cutoff must be positive"
        )

    if interval < 1:
        raise ValueError(
            "interval must be at least 1"
        )

    if workers is not None and workers < 1:
        raise ValueError(
            "workers must be at least 1"
        )

    universe = mda.Universe(psf, traj)
    all_atoms = universe.select_atoms(sel)

    if len(all_atoms) == 0:
        raise ValueError(
            f"No atoms found for selection: {sel}"
        )

    aromatic_count = sum(
        1
        for residue in all_atoms.residues
        if residue.resname in RESIDUE_TO_RING_ATOMS
    )

    cations = all_atoms.select_atoms(
        CATION_SELECTION
    )

    if aromatic_count == 0:
        raise ValueError(
            "Unable to find aromatic residues"
        )

    if aromatic_count + len(cations) < 2:
        raise ValueError(
            "The total number of aromatic residues and "
            "cation atoms is less than two"
        )

    n_frames = len(universe.trajectory)
    analyzed_frames = len(
        range(0, n_frames, interval)
    )

    print("\n============= PI STACKING INFO =============")
    print(f"Frames            : {n_frames}")
    print(f"Analyzed          : {analyzed_frames}")
    print(f"Aromatic residues : {aromatic_count}")
    print(f"Cation atoms      : {len(cations)}")
    print(f"Pi-Pi cutoff      : {pi_pi_dist_cutoff}")
    print(f"Pi-cation cutoff  : {pi_cation_dist_cutoff}")
    print(f"Interval          : {interval}")
    context = RuntimeContext.detect_desktop()
    plan = RuntimeScheduler(context).create_plan(
        task_count=n_frames or None,
        n_workers=workers,
    )
    print(f"Backend           : {plan.backend}")
    print(f"Strategy          : {plan.strategy}")
    print(f"Workers           : {plan.n_workers}")
    print("============================================\n")

    chunks = chunk_frames(
        n_frames=n_frames,
        n_workers=plan.n_workers,
    )

    tasks = [
        (
            psf,
            traj,
            sel,
            pi_pi_dist_cutoff,
            pi_cation_dist_cutoff,
            interval,
            start,
            stop,
            debug,
        )
        for start, stop in chunks
    ]

    start_time = time.perf_counter()

    partial_results = RuntimeExecutor(plan=plan).run(
        pi_stacking_worker,
        tasks,
    )

    pi_stacking = merge_events(
        partial_results
    )

    elapsed = time.perf_counter() - start_time

    print(
        f"Pi-stacking time: {elapsed:.3f} sec"
    )

    sorted_events = sorted(
        pi_stacking.items(),
        key=lambda item: (
            -len(item[1]),
            item[0],
        ),
    )

    with sta.resolve_file(
        out,
        "w",
    ) as outfile:
        header(outfile)

        for pair, frames in sorted_events:
            residue1, residue2 = pair

            frame_text = " ".join(
                str(frame)
                for frame in sorted(frames)
            )

            print(
                residue1,
                residue2,
                frame_text,
                file=outfile,
            )


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=f"stanalyzer {ANALYSIS_NAME}"
    )

    sta.add_project_args(
        parser,
        "psf",
        "traj",
        "out",
        "interval",
    )

    parser.add_argument(
        "--pi-pi-dist-cutoff",
        type=float,
        metavar="N",
        default=6.0,
        help="Distance cutoff between aromatic ring centers.",
    )

    parser.add_argument(
        "--pi-cation-dist-cutoff",
        type=float,
        metavar="N",
        default=6.0,
        help="Distance cutoff between aromatic rings and cations.",
    )

    parser.add_argument(
        "--sel",
        metavar="selection",
        default="all",
        help="Restrict the analysis to selected atoms.",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Maximum workers; defaults to automatic runtime selection.",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print per-frame event counts.",
    )

    return parser

def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(
            sta.get_settings(ANALYSIS_NAME)
        )

    write_pi_stacking(**settings)


if __name__ == "__main__":
    main()
