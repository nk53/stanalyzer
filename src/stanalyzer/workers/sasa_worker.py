import os
import tempfile
import warnings

import freesasa
import MDAnalysis as mda


def make_sasa_parameters(
    algorithm: str,
    probe_radius: float,
) -> freesasa.Parameters:
    """Create and validate FreeSASA parameters."""

    if probe_radius <= 0:
        raise ValueError("probe_radius must be positive")

    if algorithm == "shrake":
        return freesasa.Parameters(
            {
                "algorithm": freesasa.ShrakeRupley,
                "probe-radius": probe_radius,
            }
        )

    if algorithm == "lee":
        return freesasa.Parameters(
            {
                "algorithm": freesasa.LeeRichards,
                "probe-radius": probe_radius,
            }
        )

    raise ValueError(
        f"Unknown algorithm '{algorithm}'. "
        "Expected 'shrake' or 'lee'."
    )


def sasa_worker(task):
    """
    Compute SASA values for one frame chunk.

    Parameters
    ----------
    task : tuple
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

    Returns
    -------
    list[tuple[int, float]]
        Frame number and SASA value for each analyzed frame.
    """

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
    ) = task

    universe = mda.Universe(psf, traj)
    atom_group = universe.select_atoms(sel)

    if len(atom_group) == 0:
        raise ValueError(
            f"No atoms found for selection: {sel}"
        )

    parameters = make_sasa_parameters(
        algorithm=algorithm,
        probe_radius=probe_radius,
    )

    freesasa.setVerbosity(
        freesasa.normal
        if debug
        else freesasa.nowarnings
    )

    rows: list[tuple[int, float]] = []

    # Ensure interval is based on global frame indices rather than
    # restarting independently at each chunk boundary.
    first_frame = start + (-start % interval)

    for frame_idx in range(
        first_frame,
        stop,
        interval,
    ):
        ts = universe.trajectory[frame_idx]
        temp_name = None

        try:
            with tempfile.NamedTemporaryFile(
                suffix=".pdb",
                delete=False,
            ) as tmp:
                temp_name = tmp.name

            with warnings.catch_warnings():
                if not debug:
                    warnings.simplefilter(
                        "ignore",
                        UserWarning,
                    )

                atom_group.write(temp_name)

            structure = freesasa.Structure(
                temp_name
            )

            result = freesasa.calc(
                structure,
                parameters,
            )

            total_sasa = float(
                result.totalArea()
            )

            rows.append(
                (
                    int(ts.frame),
                    total_sasa,
                )
            )

            if debug:
                print(
                    f"Frame {ts.frame}: "
                    f"SASA={total_sasa:.5f} A^2"
                )

        finally:
            if temp_name is not None:
                try:
                    os.remove(temp_name)
                except FileNotFoundError:
                    pass

    return rows