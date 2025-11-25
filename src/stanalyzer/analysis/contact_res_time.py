import argparse
from typing import Optional, TypeAlias, cast

import MDAnalysis as mda
import numpy as np

import stanalyzer.cli.stanalyzer as sta
from stanalyzer.cli.validators import p_int, p_float

ANALYSIS_NAME = 'contact_res_time'
CONTACT: TypeAlias = tuple[str, int, str, int]  # RESNAME-RESID RESNAME-RESID of contact pair
FRAME_COUNTS: TypeAlias = list[int]             # number of frames for a single contact event
FRAME_DIST: TypeAlias = tuple[float, float]     # mean +/- standard deviation


def header(outfile: sta.FileLike | None = None, np_formatted: bool = False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "residue1_name residue1_ID residue2_name residue2_ID mean_time std"
    else:
        header_str = "#residue1_name residue1_ID residue2_name residue2_ID mean_time std"

    print(header_str, file=outfile)

    return header_str


def write_contact_res_time(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
                           out: sta.FileRef, threshold: float = 5.0,
                           ref_psf: Optional[sta.FileRef] = None,
                           ref_coor: Optional[sta.FileRef] = None, ref_frame_num: int = 1,
                           interval: int = 1) -> None:
    """Writes RMSD to `out` file"""

    if ref_psf is None:
        ref_psf = cast(sta.FileRef, psf)
    if ref_coor is None:
        ref_coor = cast(sta.FileRef, traj[0])

    universe = mda.Universe(psf, traj)

    # Initialize dictionaries to track contact frames and residence times
    contact_frames: dict[CONTACT, int] = {}
    residence_times: dict[CONTACT, FRAME_COUNTS] = {}
    residence_dist: dict[CONTACT, FRAME_DIST] = {}

    # Iterate over frames in the trajectory
    for step_num, ts in enumerate(universe.trajectory, start=1):
        if step_num % interval:
            continue

        all_atoms = universe.select_atoms(sel)
        print(all_atoms)
        # Group residues from selected atoms
        residues = all_atoms.residues

        # Check if we have any residues
        if len(residues) == 0:
            print("No residues found in the selection.")
            continue

        # Iterate over pairs of residues
        for i, res_i in enumerate(residues):
            for j in range(i + 1, len(residues)):
                res_j = residues[j]
                com_i = res_i.atoms.center_of_mass()
                com_j = res_j.atoms.center_of_mass()

                # Calculate distance
                dist = np.linalg.norm(com_i - com_j)
                contact = (res_i.resname, res_i.resid, res_j.resname, res_j.resid)

                if dist < threshold:
                    # Track contact frame
                    contact_frames[contact] = contact_frames.get(contact, 0) + 1

                else:
                    # Reset count if contact is lost
                    if contact in contact_frames:
                        residence_times.setdefault(contact, [])
                        residence_times[contact].append(contact_frames.pop(contact))

    # Add remaining contact counts to residence times
    for contact, frames in contact_frames.items():
        residence_times.setdefault(contact, [])
        residence_times[contact].append(frames)

    # Get residence time distribution as mean ± standard deviation
    for contact, event_len in residence_times.items():
        if event_len:
            events = np.array(event_len)
            residence_dist[contact] = events.mean(), events.std()
    else:
        residence_dist[contact] = (0.0, 0.0)  # Handle case with no events

    with sta.resolve_file(out, 'w') as outfile:
        print("# Residue1_Resname Residue1_ID Residue2_Resname Residue2_ID Mean_Time Std",
              file=outfile)
        for contact, (mean, std) in residence_dist.items():
            # res_i_name, res_i_id, res_j_name, res_j_id = contact
            print(*contact, f"{mean:.5g} {std:.5g}", file=outfile)

    print(f"Residence times written to {out}")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval')
    parser.add_argument('-rp', '--ref-psf', '--ref-psf-path', type=sta.InputFile,
                        metavar='FILE',
                        help="PSF to use for reference, if not same as --psf")
    parser.add_argument('-rc', '--ref-coor', '--ref-coor-path', type=sta.InputFile,
                        metavar='FILE',
                        help="Coordinate file to use for reference, if not same as --traj")
    parser.add_argument('--sel', metavar='selection',
                        help="Atom selection for RMSD calculation")
    parser.add_argument('--threshold', type=p_float, metavar='N', default='5.0',
                        help="Distance cutoff for calculating the residence time.")
    parser.add_argument('-rn', '--ref-frame-num', type=p_int, default=1, metavar='N',
                        help="Frame to use for reference coordinates (default: 1). "
                        "Only meaningful if --ref-frame-type is 'specific'")

    return parser


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contact_res_time(**settings)


if __name__ == '__main__':
    main()
