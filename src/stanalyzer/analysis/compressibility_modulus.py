import argparse
import typing as t

import stanalyzer.cli.stanalyzer as sta
import numpy as np
import MDAnalysis as mda

ANALYSIS_NAME = 'compressibility_modulus'

T = t.TypeVar('T')
Array1D: t.TypeAlias = list[T]
Array2D: t.TypeAlias = list[Array1D[T]]
FL: t.TypeAlias = list[float]
Tup5: t.TypeAlias = tuple[FL, FL, FL, FL, float]


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "Area compressibility_modulus from box size fluctuation"

    print(header_str, file=outfile)

    return header_str


def write_compressibility_modulus(psf: sta.FileRef, traj: sta.FileRefList,
                                  out: sta.FileRef, temp: float) -> None:
    """Writes Compressibility Modulus to 'out' file"""

    # estimate KA from fluctuation of Box
    kB = 1.380649e-23  # J/K ; Boltzmann constant
    kBT = kB*temp      # J=N*m = 10^5 dyn * 10^2 cm
    kBT = kBT*1.0e7    # in units of dyn * cm
    ang2cm = 1.e-8     # cm/A

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed - All frames
    framenum = int(u.trajectory.n_frames)

    # set up box_area
    sa = np.zeros([framenum], dtype=float)

    for i in range(0, framenum):
        if framenum < 100:
            print(f'# processing {i+1}/{framenum}')
        else:
            if i % int(framenum/100):
                print(f'# processing {i+1}/{framenum}')
        ts = u.trajectory[i]

        # get box size: Assumes an orthorhombic box
        sa[i] = ts.dimensions[0]*ts.dimensions[1]

    # Mean xy_area & its variance
    asa = np.average(sa)  # mean
    vsa = np.var(sa)      # simple variance

    # area compressibility modulus
    # KA = kBT * <A>/(<A^2>-<A>^2))
    # KA = kBT * ave_area / var_are
    ka = kBT * asa/vsa   # in units of  dyn*cm/A^2
    ka /= ang2cm*ang2cm  # in the unit of dyn/cm

    sout = f'K_A = {ka:10.5f} (dyn/cm)\n'

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        print(sout, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    parser.add_argument('--temp', type=float, help='System Temperature')

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_compressibility_modulus(**settings)


if __name__ == '__main__':
    main()
