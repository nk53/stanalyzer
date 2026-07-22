import argparse
import time
from collections import defaultdict

import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

from stanalyzer.workers.salt_bridge_worker import (
    DEFAULT_NEGATIVE_DEF,
    DEFAULT_POSITIVE_DEF,
    build_selection,
    salt_bridge_worker,
)
from stanalyzer.runtime import (
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    chunk_frames,
)


ANALYSIS_NAME = "salt_bridge"

SaltBridge = tuple[str, str]
SaltBridgeFrames = dict[SaltBridge, set[int]]


def header(
    outfile: sta.FileLike | None = None,
) -> str:
    header_str = "#residue1 residue2 frames"

    if outfile is not None:
        print(header_str, file=outfile)

    return header_str


def merge_salt_bridges(
    partial_results: list[SaltBridgeFrames],
) -> SaltBridgeFrames:
    """Union frame sets returned by all workers."""

    merged: dict[SaltBridge, set[int]] = defaultdict(set)

    for partial in partial_results:
        for bridge, frames in partial.items():
            merged[bridge].update(frames)

    return dict(merged)


def write_salt_bridge(
    psf: sta.FileRef,
    traj: sta.FileRefList,
    out: sta.FileRef,
    positive_sel: str = "all",
    negative_sel: str = "all",
    positive_def: str | None = None,
    negative_def: str | None = None,
    dist_cutoff: float = 4.5,
    interval: int = 1,
    workers: int | None = None,
    debug: bool = False,
) -> None:
    """Calculate residue-level salt bridges over a trajectory."""

    if dist_cutoff <= 0:
        raise ValueError("dist_cutoff must be positive")

    if interval < 1:
        raise ValueError("interval must be at least 1")

    if workers is not None and workers < 1:
        raise ValueError("workers must be at least 1")

    positive_selection = build_selection(
        outer_selection=positive_sel,
        atom_definition=positive_def,
        default_definition=DEFAULT_POSITIVE_DEF,
    )

    negative_selection = build_selection(
        outer_selection=negative_sel,
        atom_definition=negative_def,
        default_definition=DEFAULT_NEGATIVE_DEF,
    )

    # Parent Universe is only for validation and trajectory metadata.
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

    n_frames = len(universe.trajectory)
    analyzed_frames = len(range(0, n_frames, interval))

    print("\n============ SALT BRIDGE INFO ============")
    print(f"Frames          : {n_frames}")
    print(f"Analyzed        : {analyzed_frames}")
    print(f"Acidic atoms    : {len(acidic)}")
    print(f"Basic atoms     : {len(basic)}")
    print(f"Distance cutoff : {dist_cutoff} Å")
    print(f"Interval        : {interval}")
    context = RuntimeContext.detect()
    plan = RuntimeScheduler(context).create_plan(
        task_count=n_frames or None,
        n_workers=workers,
    )
    print(f"Backend         : {plan.backend}")
    print(f"Strategy        : {plan.strategy}")
    print(f"Workers         : {plan.n_workers}")
    print("==========================================\n")

    chunks = chunk_frames(
        n_frames=n_frames,
        n_workers=plan.n_workers,
    )

    tasks = [
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
        for start, stop in chunks
    ]

    start_time = time.perf_counter()

    partial_results = RuntimeExecutor(plan=plan).run(
        salt_bridge_worker,
        tasks,
    )

    salt_bridges = merge_salt_bridges(
        partial_results
    )

    elapsed = time.perf_counter() - start_time

    print(f"Salt bridge time: {elapsed:.3f} sec")

    # Sort by number of observed frames, then by residue pair to make
    # output deterministic across different worker counts.
    sorted_bridges = sorted(
        salt_bridges.items(),
        key=lambda item: (
            -len(item[1]),
            item[0],
        ),
    )

    with sta.resolve_file(out, "w") as outfile:
        header(outfile)

        for bridge, frames in sorted_bridges:
            residue1, residue2 = bridge
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

    print(f"Salt bridge results written to {out}")


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
        "--dist-cutoff",
        type=p_float,
        metavar="N",
        default=4.5,
        help="Distance cutoff between oppositely charged atoms.",
    )

    parser.add_argument(
        "--positive-sel",
        metavar="selection",
        default="all",
        help="Scope used when selecting positively charged atoms.",
    )

    parser.add_argument(
        "--negative-sel",
        metavar="selection",
        default="all",
        help="Scope used when selecting negatively charged atoms.",
    )

    parser.add_argument(
        "--positive-def",
        metavar="selection",
        default=None,
        help="Custom definition of positively charged atoms.",
    )

    parser.add_argument(
        "--negative-def",
        metavar="selection",
        default=None,
        help="Custom definition of negatively charged atoms.",
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
        help="Print per-frame bridge counts.",
    )

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(
            sta.get_settings(ANALYSIS_NAME)
        )

    write_salt_bridge(**settings)


if __name__ == "__main__":
    main()
