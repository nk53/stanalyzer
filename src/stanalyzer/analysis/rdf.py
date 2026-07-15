import argparse
import warnings
import time

import MDAnalysis as mda
import numpy as np

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_float

from stanalyzer.workers.rdf_worker import rdf_worker
from stanalyzer.runtime import (
    RuntimeContext,
    RuntimeExecutor,
    RuntimeScheduler,
    chunk_frames,
)


if hasattr(np, "VisibleDeprecationWarning"):
    warnings.simplefilter(
        "ignore",
        category=np.VisibleDeprecationWarning,
    )

ANALYSIS_NAME = "Radial distribution function"


def header(outfile=None):
    header_str = "#Bin RDF"
    print(header_str, file=outfile)
    return header_str


def box_volume(dimensions):
    x, y, z = dimensions[:3]

    if len(dimensions) < 6:
        return x * y * z

    alpha, beta, gamma = np.deg2rad(dimensions[3:6])
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)

    volume_factor = np.sqrt(np.clip(
        1.0
        - cos_alpha**2
        - cos_beta**2
        - cos_gamma**2
        + 2.0 * cos_alpha * cos_beta * cos_gamma,
        0.0,
        None,
    ))

    return x * y * z * volume_factor


def write_rdf(
    psf,
    traj,
    out,
    sel1,
    sel2,
    bin_size,
    step=1,
    method="bruteforce",
    workers=None,
):
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")

    if step < 1:
        raise ValueError("step must be at least 1")

    if workers is not None and workers < 1:
        raise ValueError("workers must be at least 1")

    u = mda.Universe(psf, traj)

    sel_atoms1 = u.select_atoms(sel1)
    sel_atoms2 = u.select_atoms(sel2)

    if len(sel_atoms1) == 0:
        raise ValueError(f"No atoms found for selection 1: {sel1}")

    if len(sel_atoms2) == 0:
        raise ValueError(f"No atoms found for selection 2: {sel2}")

    same_atoms = np.array_equal(
        sel_atoms1.indices,
        sel_atoms2.indices,
    )

    n_pairs = len(sel_atoms1) * len(sel_atoms2)
    if same_atoms:
        n_pairs -= len(sel_atoms1)

    if n_pairs <= 0:
        raise ValueError("RDF requires at least one atom pair")

    min_box_length = None
    volume_sum = 0.0
    analyzed_frames = 0

    n_frames = len(u.trajectory)

    for frame_idx in range(0, n_frames, step):
        ts = u.trajectory[frame_idx]

        if ts.dimensions is None:
            raise ValueError(
                "RDF requires trajectory frames with box dimensions"
            )

        frame_min = np.min(ts.dimensions[:3])
        min_box_length = (
            frame_min
            if min_box_length is None
            else min(min_box_length, frame_min)
        )
        volume_sum += box_volume(ts.dimensions)
        analyzed_frames += 1

    if analyzed_frames == 0:
        raise ValueError("No trajectory frames selected for RDF analysis")

    assert min_box_length is not None

    # Match MDAnalysis InterRDF and the v1 implementation: the largest
    # unambiguous minimum-image distance is half the shortest box length.
    rdf_range = min_box_length * 0.5

    nbins = int(rdf_range / bin_size)

    if nbins < 1:
        raise ValueError("bin_size is too large for the RDF range")

    rdf_range = nbins * bin_size

    print("\n================ RDF INFO ================")
    print(f"Frames       : {n_frames}")
    print(f"Analyzed     : {analyzed_frames}")
    print(f"Selection 1  : {len(sel_atoms1)} atoms")
    print(f"Selection 2  : {len(sel_atoms2)} atoms")
    print(f"RDF range    : {rdf_range:.3f} Å")
    print(f"Bin size     : {bin_size}")
    print(f"Bins         : {nbins}")
    print(f"Step         : {step}")
    print(f"Method       : {method}")
    context = RuntimeContext.detect_desktop()
    plan = RuntimeScheduler(context).create_plan(
        task_count=n_frames or None,
        n_workers=workers,
    )
    print(f"Backend      : {plan.backend}")
    print(f"Strategy     : {plan.strategy}")
    print(f"Workers      : {plan.n_workers}")
    print("==========================================\n")

    start_time = time.perf_counter()

    chunks = chunk_frames(
        n_frames=n_frames,
        n_workers=plan.n_workers,
    )

    tasks = [
        (
            psf,
            traj,
            sel1,
            sel2,
            rdf_range,
            nbins,
            method,
            step,
            start,
            stop,
        )
        for start, stop in chunks
    ]

    executor = RuntimeExecutor(plan=plan)

    partial_hists = executor.run(
        rdf_worker,
        tasks,
    )

    hist = np.sum(
        partial_hists,
        axis=0,
    )

    elapsed = time.perf_counter() - start_time

    print(
        f"RDF Time: {elapsed:.3f} sec"
    )

    edges = np.linspace(
        0.0,
        rdf_range,
        nbins + 1,
    )

    bins = (
        edges[:-1]
        + edges[1:]
    ) * 0.5

    shell_volumes = (
        4.0
        / 3.0
        * np.pi
        * (edges[1:] ** 3 - edges[:-1] ** 3)
    )
    average_volume = volume_sum / analyzed_frames
    density = n_pairs / average_volume
    rdf = hist / (
        analyzed_frames
        * shell_volumes
        * density
    )

    with sta.resolve_file(
        out,
        "w",
    ) as outfile:

        header(outfile)

        for r, value in zip(
            bins,
            rdf,
        ):
            print(
                f"{r:.2f} {value:.4f}",
                file=outfile,
            )


def get_parser():

    parser = argparse.ArgumentParser(
        prog=f"stanalyzer {ANALYSIS_NAME}"
    )

    sta.add_project_args(
        parser,
        "psf",
        "traj",
        "out",
    )

    parser.add_argument(
        "-sel1",
        metavar="select",
        help="Reference selection",
    )

    parser.add_argument(
        "-sel2",
        metavar="select",
        help="Target selection",
    )

    parser.add_argument(
        "-bin-size",
        type=p_float,
        default=0.1,
        help="Bin size",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Analyze every Nth frame",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Maximum workers; defaults to automatic runtime selection",
    )

    parser.add_argument(
        "--method",
        choices=[
            "bruteforce",
            "nsgrid",
            "pkdtree",
        ],
        default="bruteforce",
    )

    return parser


def main(settings=None):

    if settings is None:
        settings = dict(
            sta.get_settings(
                ANALYSIS_NAME
            )
        )

    write_rdf(**settings)


if __name__ == "__main__":
    main()
