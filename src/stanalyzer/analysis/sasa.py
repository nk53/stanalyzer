import argparse
import time

import MDAnalysis as mda

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

from stanalyzer.workers.sasa_worker import (
    make_sasa_parameters,
    sasa_worker,
)
from stanalyzer.runtime import (
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    chunk_frames,
)


ANALYSIS_NAME = "sasa"


def header(
    outfile: sta.FileLike | None = None,
) -> str:
    header_str = "# Frame SASA"

    if outfile is not None:
        print(header_str, file=outfile)

    return header_str


def write_sasa(
    psf: sta.FileRef,
    traj: sta.FileRefList,
    sel: str,
    probe_radius: float,
    algorithm: str,
    out: sta.FileRef,
    interval: int = 1,
    workers: int | None = None,
    debug: bool = False,
) -> None:
    """
    Calculate solvent-accessible surface area for selected atoms.

    Frames are divided into chunks and processed independently.
    Each worker creates its own MDAnalysis Universe and returns
    rows of:

        (frame_number, sasa_value)
    """

    if interval < 1:
        raise ValueError(
            "interval must be at least 1"
        )

    if workers is not None and workers < 1:
        raise ValueError(
            "workers must be at least 1"
        )

    # Validate algorithm and radius before starting workers.
    make_sasa_parameters(
        algorithm=algorithm,
        probe_radius=probe_radius,
    )

    # Parent Universe is used for metadata and early validation.
    universe = mda.Universe(
        psf,
        traj,
    )

    atom_group = universe.select_atoms(
        sel
    )

    if len(atom_group) == 0:
        raise ValueError(
            f"No atoms found for selection: {sel}"
        )

    n_frames = len(
        universe.trajectory
    )

    analyzed_frames = len(
        range(
            0,
            n_frames,
            interval,
        )
    )

    print("\n=============== SASA INFO ================")
    print(f"Frames       : {n_frames}")
    print(f"Analyzed     : {analyzed_frames}")
    print(f"Selection    : {sel}")
    print(f"Atoms        : {len(atom_group)}")
    print(f"Algorithm    : {algorithm}")
    print(f"Probe radius : {probe_radius}")
    print(f"Interval     : {interval}")
    context = RuntimeContext.detect_desktop()
    plan = RuntimeScheduler(context).create_plan(
        task_count=n_frames or None,
        n_workers=workers,
    )
    print(f"Backend      : {plan.backend}")
    print(f"Strategy     : {plan.strategy}")
    print(f"Workers      : {plan.n_workers}")
    print("==========================================\n")

    chunks = chunk_frames(
        n_frames=n_frames,
        n_workers=plan.n_workers,
    )

    tasks = [
        (
            psf,
            traj,
            sel,
            probe_radius,
            algorithm,
            interval,
            debug,
            start,
            stop,
        )
        for start, stop in chunks
    ]

    start_time = time.perf_counter()

    executor = RuntimeExecutor(plan=plan)

    partial_rows = executor.run(
        sasa_worker,
        tasks,
    )

    # Merge worker outputs.
    rows = [
        row
        for partial in partial_rows
        for row in partial
    ]

    # ProcessPoolExecutor results may arrive in a different order
    # if executor implementation changes later.
    rows.sort(
        key=lambda row: row[0]
    )

    elapsed = (
        time.perf_counter()
        - start_time
    )

    print(
        f"SASA Time: {elapsed:.3f} sec"
    )

    with sta.resolve_file(
        out,
        "w",
    ) as outfile:
        header(outfile)

        for frame, total_sasa in rows:
            print(
                f"{frame} {total_sasa:.5f}",
                file=outfile,
            )

    print(
        f"SASA results written to {out}"
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
        "--sel",
        metavar="selection",
        required=True,
        help="Atom selection for SASA calculation.",
    )

    parser.add_argument(
        "--probe-radius",
        type=p_float,
        default=1.4,
        help="Probe radius used by FreeSASA.",
    )

    parser.add_argument(
        "--algorithm",
        choices=[
            "shrake",
            "lee",
        ],
        default="shrake",
        help="FreeSASA algorithm.",
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
        help="Print individual frame results.",
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

    write_sasa(**settings)


if __name__ == "__main__":
    main()
