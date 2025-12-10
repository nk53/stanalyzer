import argparse
import re
import typing as t
import MDAnalysis as mda
import MDAnalysis.transformations as transformations
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from scipy.spatial import Voronoi

import stanalyzer.cli.stanalyzer as sta
from . import mol_atomgroup_util as mymol
from . import leaflet_util as myleaflet
from . import voronoi_analysis as myvorn
from .voronoi_analysis import (
    List2D,
    NDFloat64,
    NDInt64,
)

ANALYSIS_NAME = 'voronoi_shell_comp'


# LeafletAssignmentMethod: t.TypeAlias = t.Literal['mda', 'zpos']
OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]


# --- The following are hard set for membrane analysis
dim = 2              # 2 dimension
nimage = int(3**dim)  # number of total primary+images for 2D

nside = 2     # up/dn
sside = ["up", "dn"]


def process_args(sel: str, split_to_mol: str = '') -> ProcessedArgs:
    """
    ----------
    Process arguments
    ----------
    """

    # selections for individual molecule types
    selection = re.split(';|,', f'{sel:s}')
    ntype = len(selection)                  # number of unique molecule types
    for i in range(0, ntype):
        selection[i] = selection[i].strip()

    # in case split_to_mol is None
    if not split_to_mol:
        nsplit = 0
    # normal process of split_to_mol
    else:
        split = re.split(';|,', f'{split_to_mol:s}')
        nsplit = len(split)
        for i in range(0, nsplit):
            split[i] = split[i].strip()

    qsplit = []
    if nsplit < ntype:  # add more qsplit options
        for i in range(0, nsplit):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)
        for i in range(nsplit, ntype):
            qsplit.append(True)  # default values
    else:  # get split upto ntype
        qsplit = []
        for i in range(0, ntype):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)

    return ProcessedArgs(selection, ntype, qsplit)


def write_ave_std_leaflet(nside: int, sside: list[str],
                          ntype: int, name_type: list[str],
                          SLnmol_type: List2D | NDInt64,
                          nshell: int, array1: NDFloat64, array2: NDFloat64,
                          odir: str, obstype: str, suffix: str) -> None:
    # sside: leaflet name
    # array1: average
    # array2: std
    # obstype: output type name
    for i in range(0, nside):
        tnmol_type = SLnmol_type[i]  # component numbers for the leaflet
        side = sside[i]
        for j in range(0, ntype):
            tnamej = name_type[j].strip("*")

            # handle nmol_type = 0 case
            if tnmol_type[j] == 0:
                sout = f'no {tnamej} in {side} leaflet.\n'
                sta.write_to_outfile(
                    f'{odir}/NA_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat',
                    sout)
                continue

            # Write output string for nmol_type > 0
            sout = f'# {side} {tnamej}: shell [{obstype} & std]_shell ntype= {ntype}'
            sout += f' nshell= {nshell}\n'
            sout += '#     shell'
            for m in range(0, ntype):
                tnamem = name_type[m].strip("*")
                sout += f' {tnamem:21s}'
            sout += '\n'
            for k in range(1, nshell+1):
                sout += f' {k:10d}'  # shell number
                for m in range(0, ntype):
                    sout += f' {array1[i, j, k, m]:10.5f} {array2[i, j, k, m]:10.5f}'
                sout += '\n'
            sta.write_to_outfile(
                f'{odir}/ave_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat', sout)
            # print(sout)


def write_ave_std_bilayer(ntype: int, name_type: list[str],
                          Snmol_type: list[int] | NDInt64,
                          nshell: int, array1: NDFloat64, array2: NDFloat64,
                          odir: str, obstype: str, suffix: str) -> None:
    # array1: average
    # array2: std
    # obstype: output type name
    side = "bilayer"
    for j in range(0, ntype):
        tnamej = name_type[j].strip("*")

        # handle nmol_type = 0 case
        if Snmol_type[j] == 0:
            sout = f'no {tnamej} in the bilayer.\n'
            sta.write_to_outfile(
                f'{odir}/NA_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat',
                sout)
            continue

        # Write output for nmol_type > 0
        sout = f'# {side} {tnamej}: shell [{obstype} & std]_shell ntype= {ntype} '
        sout += f'nshell= {nshell}\n'
        sout += '#     shell'
        for m in range(0, ntype):
            tnamem = name_type[m].strip("*")
            sout += f' {tnamem:21s}'
        sout += '\n'
        for k in range(1, nshell+1):
            sout += f' {k:10d}'  # shell number
            for m in range(0, ntype):
                sout += f' {array1[j, k, m]:10.5f} {array2[j, k, m]:10.5f}'
            sout += '\n'
        sta.write_to_outfile(
            f'{odir}/ave_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat', sout)
        # print(sout)


def write_time_series_leaflet(framenum: int, interval: int, time_step: float,
                              nside: int, sside: list[str],
                              ntype: int, name_type: list[str],
                              SLnmol_type: List2D | NDInt64,
                              nshell: int, array: NDFloat64,
                              odir: str, obstype: str, suffix: str) -> None:
    # sside: leaflet name
    # array: time series
    # obstype: output type name
    for j in range(0, nside):
        side = sside[j]
        tnmol_type = SLnmol_type[j]

        for k in range(0, ntype):
            tnamek = name_type[k].strip("*")

            # handle nmol_type = 0 case
            if tnmol_type[k] == 0:
                sout = f'no {tnamek} in {side} leaflet.\n'
                sta.write_to_outfile(
                    f'{odir}/NA_{side}_{tnamek.lower()}_{obstype}_{suffix}.dat',
                    sout)
                continue

            # write output for nmol_type > 0
            sout = f'# {side} {tnamek}: time series of shell compositions\n'
            sout += f'# time [{obstype}]_shell ntype= {ntype} nshell= {nshell}\n'
            sout += '#          '
            for m in range(1, nshell+1):
                tsshell = f'shell{m}'
                sout += f'{tsshell:44s}'
            sout += '\n'
            sout += '#          '
            for m in range(1, nshell+1):
                for n in range(0, ntype):
                    tnamen = name_type[n].strip("*")
                    sout += f' {tnamen:10s}'
            sout += '\n'
            for i in range(0, framenum):
                sout += f' {time_step*(interval*i+1):10.5f}'
                # sout += f' {interval*i+1:10d}'
                for m in range(1, nshell+1):
                    for n in range(0, ntype):
                        sout += f' {array[i, j, k, m, n]:10.5f}'
                sout += '\n'
            fout = f'{odir}/time_{side}_{tnamek.lower()}_{obstype}_{suffix}.dat'
            sta.write_to_outfile(fout, sout)
            # print(sout)


def write_time_series_bilayer(framenum: int, interval: int, time_step: float,
                              ntype: int, name_type: list[str],
                              Snmol_type: list[int] | NDInt64,
                              nshell: int, array: NDFloat64,
                              odir: str, obstype: str, suffix: str) -> None:
    # array: time series
    # obstype: output type name
    side = "bilayer"
    for j in range(0, ntype):
        tnamej = name_type[j].strip("*")

        # handle nmol_type = 0 case
        if Snmol_type[j] == 0:
            sout = f'no {tnamej} in the bilayer.\n'
            sta.write_to_outfile(
                f'{odir}/NA_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat',
                sout)
            continue

        # write output for nmol_type > 0
        sout = f'# {side} {tnamej}: time series of shell compositions\n'
        sout += f'# time [{obstype}_comp]_shell ntype= {ntype} nshell= {nshell}\n'
        sout += '#          '
        for k in range(1, nshell+1):
            tsshell = f'shell{k}'
            sout += f'{tsshell:44s}'
        sout += '\n'
        sout += '#          '
        for k in range(1, nshell+1):
            for m in range(0, ntype):
                tnamem = name_type[m].strip("*")
                sout += f' {tnamem:10s}'
        sout += '\n'
        for i in range(0, framenum):
            sout += f' {time_step*(interval*i+1):10.5f}'
            # sout += f' {interval*i+1:10d}'
            for k in range(1, nshell+1):
                for m in range(0, ntype):
                    sout += f' {array[i, j, k, m]:10.5f}'
            sout += '\n'
        fout = f'{odir}/time_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat'
        sta.write_to_outfile(fout, sout)
        # print(sout)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval', 'time_step')
    parser.add_argument('--sel', metavar='selection', default='',
                        help='Selection for individual molecule type with a format, '
                        'segid/resname/moltype MOLECULENAME and name ATOMNAMES. For the CHARMM '
                        'format topology, an selection for POPC can be "resname POPC and (name C2 '
                        'or name C21 or name C31)".')
    parser.add_argument('--split', action='store', default='',
                        help='Y/N. If N, the atom group for the selection is considered as that '
                        'for a single molecule. If Y, the atom group is further splitted to '
                        'molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--sel-sys', metavar='selection', default='',
                        help='Selection for system atom groups in membranes for bilayer '
                        'recentering.')
    # ## Commented out: Activate when non-plana bilayers are supported
    # parser.add_argument('--lam', metavar='OPT', default='mda', choices=['mda', 'zpos'],
    #                     help='Leleat assignment method. mda: MDAnalysis.analysis.leaflet; '
    #                     'zpos: Z-position. Default: mda')
    parser.add_argument('--suffix', type=str, default='0',
                        help='Suffix  to output file(s)')
    parser.add_argument('--otype', metavar='OPT', default='outl', choices=['outl', 'outb'],
                        help="Output type. outl: leaflets; outb: bilayer. Default: mda")
    parser.add_argument('--qt', action='store_true', default=False,
                        help='Write output time series if provided')
    parser.add_argument('--qa', action='store_true', default=False,
                        help='Write output averages if provided')
    return parser


def run_voronoi_shell_comp(sel: str, split: str, sel_sys: str, qt: bool, qa: bool,
                           psf: sta.FileRef, traj: sta.FileRefList, time_step: float | str,
                           suffix: str,
                           interval: int = 1,  # lam: LeafletAssignmentMethod='mda',
                           otype: OutputFileType = 'outl') -> None:
    """
    ----------
    Calculate shell compositions around inidividual molecule types.

    ----------
    """

    # process non-system arguments
    selection, ntype, qsplit = process_args(sel, split)
    # method=lam
    method: t.Literal['zpos', 'mda'] = 'zpos'  # use this option for planar bilayers
    outtype = otype
    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    # print summary of arguments
    for i in range(0, ntype):
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    print('Suffix to output files', suffix)
    print(f'Bilayer is recentered at z = 0 using {sel_sys}:')
    if method == 'zpos':
        print('Leaflets are assigned based on z-positions')
    elif method == 'mda':
        print('Leaflets are assgined using a modified version of MDA LeafletFinder')
    print('Output type', outtype)
    print('Write averge shell-composition numbers/fractions', qa)
    print('Write time sereis output', qt)
    print(f'Analysis will be done every {interval} frames')

    # Set output path
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
    ag_cent = u.select_atoms(sel_sys)
    ag_all = u.atoms

    workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point=origin),
                transformations.center_in_box(ag_cent, point=origin),
                transformations.unwrap(ag_all)]

    u.trajectory.add_transformations(*workflow)

    # generate molecule atom groups
    name_type0, nmol_type0, nmol0, id_type0, ag0 = \
        mymol.generate_mol_groups(u, ntype, selection, qsplit)

    # atom group generation
    atomgroup = u.atoms[[]]
    for i in range(0, nmol0):
        atomgroup += ag0[i]

    # Assumes that all lipid types are included in selections
    print('# Estimate of max_shell: from min(boxx,boxy)/2/sqrt(APL)')
    hard_max_shell = 5  # do not go beyond this shell
    tss = u.trajectory[0]  # first frame
    tbox = tss.dimensions[:2]  # x,y dimension
    tapl = 2.0*tbox[0]*tbox[1]/float(nmol0)  # APL
    tlen = np.sqrt(tapl)  # an estimate of x/y-dimension of a component
    tboxm = np.min(tbox)
    max_shell = int(tboxm/tlen/2.0)  # the esimate
    print(
        f'# estimated max_shell = {max_shell} & hard_set_max_shell = {hard_max_shell}')
    if max_shell > hard_max_shell:
        max_shell = hard_max_shell
    print(f'# max_shell: set to {max_shell}')

    # Setup time series data
    # Shell info.: list array of time series
    # numb/frac_comp[time (,leaflet), type, shell, type]...:
    numb_comp, frac_comp =\
        myvorn.setup_shell_comp_arrays(
            framenum, nside, ntype, max_shell, outtype)

    # accumuated Snmol_type & SLnmol_type for later use # Check for Empty data
    if outtype == 'outb':
        Snmol_type = np.zeros([ntype], dtype=int)
    elif outtype == 'outl':
        SLnmol_type = np.zeros([nside, ntype], dtype=int)

    # START: LOOP over frames
    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]

        # do frame-wise bilayer recentering - remaining translation
        Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
        zref = np.zeros([2], dtype=float)
        for iside in range(0, nside):
            zref[iside] = np.mean(Lag_ref[iside].positions[:, 2])
        # translation for z-centering
        tran = 0, 0, -np.mean(zref)
        ts = transformations.translate(tran)(ts)
        ts = transformations.unwrap(ag_all)(ts)

        # get box size
        xtla, xtlb, xtlc = ts.dimensions[:3]
        box = np.array([xtla, xtlb], dtype=float)

        # initialize number of individual molecule types in the bilayer
        if outtype == 'outb':
            nmol_typeb = np.zeros([ntype], dtype=float)

        # Assign leaflet & generate mol groups for membrane
        name_type, Lnmol_type, Lnmol, Lid_type, Lag = \
            mymol.generate_mol_groups_memb(
                u, nside, ntype, selection, qsplit, sside, method)

        # START: LOOP over leaflets
        for iside in range(0, nside):
            print(f'# leaflet {iside}')

            # number of molecules in the leaflet
            nmol = Lnmol[iside]
            # nmuber of molecules for individual types
            nmol_type = Lnmol_type[iside]
            # molecule type index for individual molecules
            id_type = Lid_type[iside]

            # generate VC list for individual molecules &
            # assign molecule index to individual VCs
            mol_ivc, id_mol_vc =\
                myvorn.connect_mol_and_vc(u, nmol, nmol_type, Lag[iside])

            # PRINT BASIC INFO.
            sinfo = '# SYS. INFO.\n'
            sinfo += f'\t numb. of user defined mol. types= {ntype}\t'
            for itype in range(0, ntype):
                tname = name_type[itype].strip("*")
                sinfo += f' {tname}'
            sinfo += '\n\t\t\t\t\t\t'
            for itype in range(0, ntype):
                sinfo += f' {nmol_type[itype]:4d}'
            sinfo += f' ntot= {nmol}\n'
            sinfo += f'\t system size\t\t\t\t xtla= {xtla:10.5f} xtlb= {xtlb:10.5f}'
            print(sinfo)

            # PREPARE VORONOI CENTER
            print('# VORONOI TESSELLATION')
            nprim, nvorn, vc = \
                myvorn.prepare_voronoi_centers(Lag[iside], box, dim, nimage)

            # VORONOI TESSEELATION
            vorn = Voronoi(vc)

            # GENERATE Voronoi Center CONTACT INFO.
            print('# GENERATE VC CONTACT INFO.')
            contactVC, borderVC, distVC, weightVC =\
                myvorn.generate_vc_contact_info(vorn, nprim, nimage, box)

            # GENERATE MOLECULAR CONTACT INFO.
            print('# GENERATE MOL CONTACT INFO.')
            contactMol = \
                myvorn.generate_mol_contact_info(
                    mol_ivc, contactVC, id_mol_vc, nprim)

            # SHELL-WISE NUMBER & FRACTIONAL COMPOSITIONS.
            print('# GENERATE SHELL-WISE NUMBER & FRACTIONAL COMPOSITIONS')
            for imol in range(0, nmol):
                if nmol >= 1000:
                    if imol % int(nmol/10) == 0:
                        print(f'shell analysis for molecule {imol}/{nmol}')
                itype = id_type[imol]    # molecule type index

                # calculate tnumb/tfrac_comp[max_shell,ntype] for each molecule
                # NDarrays are returned
                tnumb_comp, tfrac_comp = \
                    myvorn.calculate_shell_comp(
                        imol, nmol, contactMol, ntype, mol_ivc, id_type, max_shell)

                # update time series (un-normalized) for each molecule type
                # Accumulating numbers.
                if outtype == 'outb':
                    numb_comp[i][itype], frac_comp[i][itype] = \
                        myvorn.update_shell_comp(ntype, max_shell,
                                                 numb_comp[i][itype], frac_comp[i][itype],
                                                 tnumb_comp, tfrac_comp)
                elif outtype == 'outl':
                    numb_comp[i][iside][itype], frac_comp[i][iside][itype] = \
                        myvorn.update_shell_comp(ntype, max_shell,
                                                 numb_comp[i][iside][itype],
                                                 frac_comp[i][iside][itype],
                                                 tnumb_comp, tfrac_comp)

            if outtype == 'outb':
                # update nmol_typeb for normalization
                nmol_typeb += nmol_type
                Snmol_type += nmol_type  # update total sum
            elif outtype == 'outl':
                SLnmol_type[iside] += nmol_type  # update total sum
                numb_comp[i][iside], frac_comp[i][iside] = \
                    myvorn.normalize_raw_comp(
                        numb_comp[i][iside], frac_comp[i][iside],
                        t.cast(list[float], nmol_type))
        # normalize bilayer data
        if outtype == 'outb':
            numb_comp[i], frac_comp[i] = \
                myvorn.normalize_raw_comp(
                    numb_comp[i], frac_comp[i], nmol_typeb)
    # END: LOOP over frames

    # Process to get statistics...
    if qa:
        anumb_comp = np.average(numb_comp, axis=0)
        snumb_comp = np.std(numb_comp, axis=0)

        afrac_comp = np.average(frac_comp, axis=0)
        sfrac_comp = np.std(frac_comp, axis=0)

        print('# Write average output')
        # Do not write output for molecule with numb_type = 0
        if outtype == 'outb':
            # write NUMB_COMP outputs # 1 to max_shell
            obstype = "ncomf"
            write_ave_std_bilayer(ntype, name_type, Snmol_type, max_shell,
                                  anumb_comp, snumb_comp, odir, obstype, suffix)

            # write FRAC_COMP outputs: shell 1 to max_shell
            obstype = "fcomp"
            write_ave_std_bilayer(ntype, name_type, Snmol_type, max_shell,
                                  afrac_comp, sfrac_comp, odir, obstype, suffix)
        elif outtype == 'outl':
            # write NUMB_COMP outputs # 1 to max_shell
            obstype = "ncomp"
            write_ave_std_leaflet(nside, sside, ntype, name_type,
                                  SLnmol_type, max_shell,
                                  anumb_comp, snumb_comp, odir, obstype, suffix)

            # write FRAC_COMP outputs: shell 1 to max_shell
            obstype = "fcomp"
            write_ave_std_leaflet(nside, sside, ntype, name_type,
                                  SLnmol_type, max_shell,
                                  afrac_comp, sfrac_comp, odir, obstype, suffix)

    # Time series output
    if qt:  # pass
        print('# Write time series output')
        # dt=1.0/float(framenum) # increment in time
        if outtype == 'outb':
            obstype = "ncomp"
            write_time_series_bilayer(framenum, interval, time_step,
                                      ntype, name_type, Snmol_type,
                                      max_shell, numb_comp, odir, obstype, suffix)

            # fractional composition/shell
            obstype = "fcomp"
            write_time_series_bilayer(framenum, interval, time_step,
                                      ntype, name_type, Snmol_type,
                                      max_shell, frac_comp, odir, obstype, suffix)
        elif outtype == 'outl':
            obstype = "ncomp"
            write_time_series_leaflet(framenum, interval, time_step,
                                      nside, sside, ntype, name_type, SLnmol_type,
                                      max_shell, numb_comp, odir, obstype, suffix)

            # fractional composition/shell
            obstype = "fcomp"
            write_time_series_leaflet(framenum, interval, time_step,
                                      nside, sside, ntype, name_type, SLnmol_type,
                                      max_shell, frac_comp, odir, obstype, suffix)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_voronoi_shell_comp(**settings)


if __name__ == '__main__':
    main()
