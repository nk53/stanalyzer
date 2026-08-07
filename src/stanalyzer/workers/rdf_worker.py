import MDAnalysis as mda
import numpy as np

from MDAnalysis.lib.distances import capped_distance


def rdf_worker(task):
    """
    Compute a partial RDF histogram for a frame chunk.

    Parameters
    ----------
    task : tuple
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

    Returns
    -------
    np.ndarray
        Partial RDF histogram.
    """

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
    ) = task

    # Each process gets its own Universe
    u = mda.Universe(psf, traj)

    ag1 = u.select_atoms(sel1)
    ag2 = u.select_atoms(sel2)
    same_atoms = np.array_equal(
        ag1.indices,
        ag2.indices,
    )

    hist = np.zeros(
        nbins,
        dtype=np.int64,
    )

    bin_width = rdf_range / nbins

    for frame_idx in range(start, stop):
        if frame_idx % step:
            continue

        ts = u.trajectory[frame_idx]

        _, dist = capped_distance(
            ag1.positions,
            ag2.positions,
            max_cutoff=rdf_range,
            box=ts.dimensions,
            method=method,
        )

        if same_atoms:
            dist = dist[dist > 0.0]

        # Faster than np.histogram
        idx = (
            dist / bin_width
        ).astype(np.int32)

        valid = (
            (idx >= 0)
            &
            (idx < nbins)
        )

        hist += np.bincount(
            idx[valid],
            minlength=nbins,
        )

    return hist
