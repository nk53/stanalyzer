import argparse
import re
import typing as t
from collections.abc import Sequence

from typing import Optional
import numpy as np

import MDAnalysis as mda
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup

import stanalyzer.cli.stanalyzer as sta
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol

if t.TYPE_CHECKING:
    import numpy.typing as npt

ANALYSIS_NAME = 'density_z'

NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
NDInt64: t.TypeAlias = 'npt.NDArray[np.int64]'
OutputFileType: t.TypeAlias = t.Literal['outs', 'outc']

# --- The following are hard set for membrane analysis
nside = 2     # up/dn
sside = ["up", "dn"]


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    names: list[str]
    ntype: int


class ProcessedArgSys(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]
    sel_type: list[str]
    name_type: list[str]


def process_args(sel: str, sel_name: str) -> ProcessedArgs:
    """
    ----------
    Process arguments
    ----------
    """

    selection = re.split(';|,', f'{sel:s}')
    names = re.split(';|,', f'{sel_name:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()
        names[i] = names[i].strip()

    return ProcessedArgs(selection, names, ntype)


def process_arg_sys(sel: str) -> ProcessedArgSys:
    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)
    for i in range(0, ntype):
        selection[i] = selection[i].strip()

    qsplit = []
    for i in range(0, ntype):
        qsplit.append(True)

    # Process selection strings to extract information
    sel_type: list[str] = []        # selection type
    name_type: list[str] = []       # name of molecule type
    for i in range(0, ntype):
        tmps = selection[i].split()
        sel_type.append(tmps[0])   # segid/resname/moleculetype
        name_type.append(tmps[1])  # PROA/PRO*/DSPC/...

    return ProcessedArgSys(selection, ntype, qsplit, sel_type, name_type)


# hist =generate_histogram(ag_targ,ntarg,numb_targ, xtlc,nbin)
def generate_histogram(ag_targ: list['AtomGroup'], ntarg: int,
                       xtlc: float, nbin: int) -> NDFloat64:

    hist = np.zeros([ntarg, nbin], dtype=float)  # histogram arrray

    zmin, zmax = -0.5*xtlc, 0.5*xtlc  # noqa: F841
    dz = xtlc/float(nbin)

    # loop over target
    for i in range(0, ntarg):
        pos = ag_targ[i].positions  # target positions
        natom = len(pos)
        for j in range(0, natom):
            tz = pos[j, 2]
            inx = int((tz-zmin)/dz)
            if inx < 0:
                inx += nbin
            if inx >= nbin:
                inx -= nbin
            hist[i, inx] += 1.0
    return hist


def write_individual_output(ntarg: int, array1: NDFloat64, array2: NDFloat64,
                            axtlc: float, nbin: int, numb_targ: Sequence[int] | NDInt64,
                            odir: str, names: list[str], suffix: str) -> None:

    bw = axtlc / float(nbin)

    for i in range(0, ntarg):
        tname = names[i].strip('*')
        tnum_targ = numb_targ[i]
        if tnum_targ == 0:
            sout = f'no {tname} in the bilayer.\n'
            fout = f'{odir}/NA_{tname}_nb{nbin}_{suffix}.dat'
        else:
            tarr1 = array1[i]
            tarr2 = array2[i]

            sout = f'# Average Z-dimesion = {axtlc} ; nbin = {nbin}\n'
            sout += '# Z-position     Average density + STD\n'
            sout += f'# {"":8s} {tname:21s}\n'

            for j in range(0, nbin):
                cz = -0.5 * axtlc + bw * (j + 0.5)  # bin center
                sout += f' {cz:10.5f} {tarr1[j]:10.5f} {tarr2[j]:10.5f}\n'

            fout = f'{odir}/{tname.lower()}_nb{nbin}_{suffix}.dat'
        sta.write_to_outfile(fout, sout)


def write_combined_output(ntarg: int, array1: NDFloat64, array2: NDFloat64,
                          axtlc: float, nbin: int, numb_targ: Sequence[int] | NDInt64,
                          odir: str, names: list[str], suffix: str) -> None:

    bw = axtlc / float(nbin)

    sout = f'# Average Z-dimesion = {axtlc} ; nbin = {nbin}\n'
    sout += '# Z-position      Average density + STD\n'
    sout += f'# {"":8s}'
    for i in range(0, ntarg):
        sout += f' {names[i].strip("*"):21s}'
    sout += '/n'

    # loop over bins
    for i in range(0, nbin):
        cz = -0.5 * axtlc + bw * (i + 0.5)  # bin center
        tarr1 = array1[:, i]  # Slice at i-th bin
        tarr2 = array2[:, i]

        sout += f' {cz:10.5f}'
        for j in range(0, ntarg):
            sout += f' {tarr1[j]:10.5f} {tarr2[j]:10.5f}'
        sout += '\n'

    if ntarg == 1:
        if numb_targ[0] == 0:
            sout = f'no {names[0]} is in the bilayer.\n'
            fout = f'{odir}/NA_{names[0].lower().strip("*")}_nb{nbin}_{suffix}.dat'
        else:
            fout = f'{odir}/{names[0].lower().strip("*")}_nb{nbin}_{suffix}.dat'
    else:
        fout = f'{odir}/combined_nb{nbin}_{suffix}.dat'

    sta.write_to_outfile(fout, sout)


def get_parser() -> argparse.ArgumentParser:
    program_name = f'stanalyzer {ANALYSIS_NAME}'
    description = "Computes a number density profile for selected atom groups along the Z axis."
    parser = argparse.ArgumentParser(
        prog=program_name, description=description)
    sta.add_project_args(parser, 'psf', 'traj', 'interval')  # , 'time_step')
    parser.add_argument('--sel', metavar='selection', required=True,
                        help='Atoms group(s) to compute density along Z axis')
    parser.add_argument('--sel-name', type=str, required=True,
                        help='Names of selected atom group(s)')
    parser.add_argument('--sel-sys', metavar='selection', required=True,
                        help='Selection for system atom groups for centering')
    # ## Commented out: Activate when non-planar bilayers are supported.
    # parser.add_argument('--lam', metavar='OPT', default='mda', choices=['mda', 'zpos'],
    #                     help='Leleat assignment method. '
    #                     'mda: MDAnalysis.analysis.leaflet; zpos: Z-position. Default: mda')
    parser.add_argument('--qcent', default=False, action='store_true',
                        help='If set True, bilayer midplane will be '
                        'centered using leaflet z-density profiles.')
    parser.add_argument('--suffix', type=str, default='0',
                        help='Suffix  to output file(s)')
    parser.add_argument('--otype', metavar='OPT', default='outs', choices=['outs', 'outc'],
                        help='Output type. outs: individual; outs: combined.'
                             ' Default: individual')
    parser.add_argument('--nbin', type=int, default=100,
                        help="Number density bin number: Default 100")

    return parser


def run_density_z(sel: str, sel_name: str, sel_sys: str, qcent: bool,
                  psf: sta.FileRef, traj: sta.FileRefList, suffix: str,
                  interval: int = 1,  # time_step: float|str,
                  # lam: LeafletAssignmentMethod='mda',
                  otype: OutputFileType = 'outs', nbin: int | str = 100) -> None:

    if otype not in ('outs', 'outc'):
        raise ValueError(f"Invalid otype: '{otype}'. Expected one of: 'outs', 'outc'.")

    if isinstance(nbin, str):  # converting bin number into integer if read as string
        nbin = int(nbin.split()[0])   # default value is 100

    # if isinstance(time_step, str):
    #     time_step = float(time_step.split()[0])

    # process arguments
    selection, names, ntarg = process_args(sel, sel_name)
    # process sys arguments - To generante full atom groups for individual molecules
    selection_sys, ntype_sys, qsplit_sys, sel_type_sys, name_type_sys \
        = process_arg_sys(sel_sys)

    # method=lam
    method = 'zpos'  # use this option for plana bilayers
    outtype = otype

    # print summary of arguments
    for i in range(0, ntarg):
        print(f'#Analysis will be done for {names[i]}')
    print('Suffix to output files', suffix)
    print(f'Bilayer is recentered at z = 0 using {sel_sys}')
    if method == 'zpos':
        print('Leaflets are assigned based on z-positions')
    elif method == 'mda':
        print('Leaflets are assgined using a modified version of MDA LeafletFinder')
    print(f'Output type: {outtype}')
    print(f'Number of bins: {nbin}')
    print(f'Analysis will be done every {interval} frames')

    # output dir
    odir = "./"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    # bilayer recentering - should be done before any assignments
    # - center in the box (an atom)
    # - center in the box (atom group for system)
    # - unwrap to get connectd molecules
    origin = 0, 0, 0  # np.zeros([3],dtype=float) ; it did not work
    # ag_cent = u.select_atoms(sel_sys)
    ag_cent = u.atoms[[]]
    for itype in range(0, ntype_sys):
        ag_cent += u.select_atoms(selection_sys[itype])
    ag_all = u.atoms

    workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point=origin),
                transformations.center_in_box(ag_cent, point=origin),
                transformations.unwrap(ag_all)]

    u.trajectory.add_transformations(*workflow)

    # print(f'len(u.trajectory) = {len(u.trajectory)}')
    # print(f'framenum = {framenum}')
    # sys.exit(0)

    # Generate atom groups for density profile calculations
    ag_targ: list['AtomGroup'] = []
    numb_targ = np.zeros([ntarg], dtype=int)
    for j in range(0, ntarg):
        tag = u.select_atoms(selection[j])
        ag_targ.append(tag)
        numb_targ[j] = len(tag)

    # Define time series of density profiles for easier statistics
    density_z = np.zeros([framenum, ntarg, nbin], dtype=float)
    xtlc = np.zeros([framenum], dtype=float)  # Z-dimension of box

    # Generate ref groups for leaflet assignemnt
    if method == "zpos":
        ag_leaflet = myleaflet.assign_leaflet_zpos(u, ag_cent)
    elif method == "mda":
        ag_leaflet = myleaflet.assign_leaflet(u, ag_cent)  # noqa: F841

    if qcent:
        print('### Generation of full atom groups for membrane molecules: START')
        # Get numbers of molecules in individual lipid types in ag_cent
        # Get total number of molecules, lipid type indices, & full atom groups
        # for z-density profile calculations of two leaflets
        #
        nmol_type_sys, nmol_sys, id_type_sys, ag_full_sys = \
            mymol.generate_full_mol_groups(
                u, ntype_sys, sel_type_sys, name_type_sys, qsplit_sys)

        print('### Generation of full atom groups for membrane molecules: DONE')

    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]

        # do frame-wise bilayer recentering
        if method == 'zpos':
            Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
        elif method == 'mda':
            Lag_ref = myleaflet.assign_leaflet(u, ag_cent)
        zref = np.zeros([2], dtype=float)
        for iside in range(0, nside):
            zref[iside] = np.mean(Lag_ref[iside].positions[:, 2])
        # translation for z-centering
        tran = 0, 0, -np.mean(zref)
        ts = transformations.translate(tran)(ts)
        ts = transformations.unwrap(ag_all)(ts)

        if qcent:  # Fine-tune bilayer midplane centering
            # Assign molecules to leaflet
            ag_full_sys_leaflet = \
                mymol.assign_full_ag_leaflet_from_ref_leaflet(
                    u, ag_full_sys, Lag_ref)

            # Calculate mid-plane z-position
            zcent = myleaflet.calc_bilayer_midplane(ag_full_sys_leaflet, nside)
            # translation for z-centering
            tran = 0, 0, -zcent
            ts = transformations.translate(tran)(ts)
            ts = transformations.unwrap(ag_all)(ts)

        # get box size
        boxx, boxy, boxz = ts.dimensions[:3]
        print(f'boxz= {boxz}')
        # box = np.array([boxx, boxy, boxz], dtype=float)
        xtlc[i] = boxz

        # generate Z-histograms
        hist = generate_histogram(ag_targ, ntarg, boxz, nbin)
        # divide by slab volume -> Number density (in units of vol^-1)
        hist /= (boxx*boxy*boxz/float(nbin))
        np.copyto(density_z[i], hist)

    # get averages & standard deviations
    axtlc = np.average(xtlc)
    adensity_z = np.average(density_z, axis=0)
    sdensity_z = np.std(density_z, axis=0)

    # Write outputs - Only ave & std.
    if otype == 'outs':  # Write outputs for each target
        write_individual_output(ntarg, adensity_z, sdensity_z, axtlc,
                                nbin, numb_targ, odir, names, suffix)
    elif otype == 'outc':  # Write combined output
        write_combined_output(ntarg, adensity_z, sdensity_z, axtlc,
                              nbin, numb_targ, odir, names, suffix)


def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginning of this function
    run_density_z(**settings)


if __name__ == '__main__':
    main()
