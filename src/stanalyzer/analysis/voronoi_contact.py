import argparse
import re
import typing as t

import MDAnalysis as mda
import MDAnalysis.transformations as transformations
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from scipy.spatial import Voronoi

import stanalyzer.cli.stanalyzer as sta
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol
from . import voronoi_analysis as myvorn

ANALYSIS_NAME = 'voronoi_contact'

if t.TYPE_CHECKING:
    import numpy.typing as npt

# LeafletAssignmentMethod: t.TypeAlias = t.Literal['mda', 'zpos']
OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']

NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
NDInt64: t.TypeAlias = 'npt.NDArray[np.int64]'


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
        for i in range(0, ntype):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)

    return ProcessedArgs(selection, ntype, qsplit)


def write_ave_std_leaflet(nside: int, sside: list[str],
                          ntype: int, name_type: list[str],
                          array1: NDFloat64, array2: NDFloat64,
                          odir: str, obstype: str, suffix: str) -> None:
    # sside:  leaflet name
    # array1: average
    # array2: std
    for i in range(0, nside):
        side = sside[i]
        sout = f'# {side}: {obstype} of contacting lipid types & std\n#'
        for j in range(0, ntype):
            tname = name_type[j].strip("*")
            sout += f' {tname:21s}'
        sout += '\n'
        for j in range(0, ntype):
            for k in range(0, ntype):
                sout += f' {array1[i, j, k]:10.5f} {array2[i, j, k]:10.5f}'
            sout += '\n'
        sta.write_to_outfile(f'{odir}/ave_{side}_{obstype}_{suffix}.dat', sout)
        # print(sout)


def write_ave_std_bilayer(ntype: int, name_type: list[str],
                          array1: NDFloat64, array2: NDFloat64,
                          odir: str, obstype: str, suffix: str) -> None:
    # array1: average
    # array2: std
    side = "bilayer"
    sout = f'# {side}: {obstype} of contacting lipid types & std\n#'
    for j in range(0, ntype):
        tname = name_type[j].strip("*")
        sout += f' {tname:21s}'
    sout += '\n'
    for j in range(0, ntype):
        for k in range(0, ntype):
            sout += f' {array1[j, k]:10.5f} {array2[j, k]:10.5f}'
        sout += '\n'

    sta.write_to_outfile(f'{odir}/ave_{side}_{obstype}_{suffix}.dat', sout)
    # print(sout)


def write_time_series_leaflet(framenum: int, interval: int, time_step: float,
                              nside: int, sside: list[str],
                              ntype: int, name_type: list[str],
                              array: NDFloat64, SLnmol_type: NDInt64, odir: str,
                              obstype: str, suffix: str) -> None:
    # array: time series
    for j in range(0, nside):
        side = sside[j]
        tnmol_type = SLnmol_type[j]
        for k in range(0, ntype):
            tnamek = name_type[k].strip("*")

            # handle nmol_type = 0 case
            if tnmol_type[k] == 0:
                sout = f'no {tnamek} in {side} leaflet.\n'
                sta.write_to_outfile(
                    f'{odir}/NA_time_{side}_{tnamek.lower()}_{obstype}_comp_{suffix}.dat',
                    sout)
                continue

            # write output for nmol_type > 0
            sout = f'# {side}: {obstype} of contacting lipid types for {tnamek}\n'
            sout += '#     frame'
            for m in range(0, ntype):
                tnamem = name_type[m].strip("*")
                sout += f' {tnamem:10s}'
            sout += '\n'
            for i in range(0, framenum):
                # ct = (dt*(i+1)+(cnt-1))
                # sout += f' {ct}'
                sout += f' {time_step*(interval*i+1):10.5f}'
                # sout += f' {interval*i+1:10d}'
                for m in range(0, ntype):
                    sout += f' {array[i, j, k, m]:10.5f}'
                sout += '\n'

            sta.write_to_outfile(
                f'{odir}/time_{side}_{tnamek.lower()}_{obstype}_{suffix}.dat', sout)
            # print(sout)


def write_time_series_bilayer(framenum: int, interval: int, time_step: float,
                              ntype: int, name_type: list[str],
                              array: NDFloat64,
                              Snmol_type: NDFloat64,
                              odir: str, obstype: str, suffix: str) -> None:
    # array: time series
    side = "bilayer"
    for j in range(0, ntype):
        tnamej = name_type[j].strip("*")

        # handle nmol_type = 0 case
        if Snmol_type[j] == 0:
            sout = f'no {tnamej} in {side} leaflet.\n'
            sta.write_to_outfile(
                f'{odir}/NA_time_{side}_{tnamej.lower()}_{obstype}_comp_{suffix}.dat',
                sout)
            continue

        # write output for nmol_type > 0
        sout = f'# {side}: {obstype} of contacting lipid types for {tnamej}\n'
        sout += '#     frame'
        for k in range(0, ntype):
            tnamek = name_type[k].strip("*")
            sout += f' {tnamek:10s}'
        sout += '\n'
        for i in range(0, framenum):
            sout += f' {time_step*(interval*i+1):10.5f}'
            # sout += f' {interval*i+1:10d}'
            for k in range(0, ntype):
                sout += f' {array[i, j, k]:10.5f}'
            sout += '\n'

        sta.write_to_outfile(
            f'{odir}/time_{side}_{tnamej.lower()}_{obstype}_{suffix}.dat', sout)
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
                        help='Selection for system atom groups for membrane recentering ')
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


def run_voronoi_contact(sel: str, split: str, sel_sys: str, qt: bool, qa: bool,
                        psf: sta.FileRef, traj: sta.FileRefList, time_step: float | str,
                        suffix: str,
                        interval: int = 1,  # lam: LeafletAssignmentMethod='mda',
                        otype: OutputFileType = 'outl') -> None:
    """
    ----------
    Analyze contacts between inidividual molecul types.

    Contacts are calculated between individual molecules in individual leaflets.
    Then, avaraged for contact between individual lipid types.
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
    print('Write averge contact numbers/fractions', qa)
    print('Write time sereis output', qt)
    print('Analysis will be done every {interval} frames')

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
    ag_cent = u.select_atoms(sel_sys)
    ag_all = u.atoms

    workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point=origin),
                transformations.center_in_box(ag_cent, point=origin),
                transformations.unwrap(ag_all)]

    u.trajectory.add_transformations(*workflow)

    # generate molecule atom groups
    name_type0, nmol_type0, nmol0, id_type0, ag0 = \
        mymol.generate_mol_groups(u, ntype, selection, qsplit)

    # Atomgroup for all atoms
    atomgroup = u.atoms[[]]
    for i in range(0, nmol0):
        atomgroup += ag0[i]

    # Setup time series data
    if outtype == 'outb':  # bilayer output
        ncontact = np.zeros([framenum, ntype, ntype], dtype=float)
        fcontact = np.zeros([framenum, ntype, ntype], dtype=float)
        # normalization factor
        tnmol_type = np.zeros([framenum, ntype], dtype=float)
        # accumulated numbers for each type: Snmol_type
        Snmol_type = np.zeros([ntype], dtype=int)
    elif outtype == 'outl':  # leaflet output
        ncontact = np.zeros([framenum, nside, ntype, ntype], dtype=float)
        fcontact = np.zeros([framenum, nside, ntype, ntype], dtype=float)
        # accumuated numbers for each type: SLnmol_type
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

            # update accumulated nmol_type
            if outtype == 'outb':
                Snmol_type += nmol_type  # update total sum
            elif outtype == 'outl':
                SLnmol_type[iside] += nmol_type  # update total sum

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

            # CALCULATE CONTACTS BTWEEN MOLECULE TYPES  (average number of contacting lipid types)
            tncontact, tfcontact =\
                myvorn.calculate_contact(
                    nmol, ntype, nmol_type, mol_ivc, id_type, contactMol)

            # update time series data
            if outtype == 'outb':
                ncontact[i] += (tncontact.T * nmol_type).T
                tnmol_type[i] += nmol_type
            elif outtype == 'outl':
                np.copyto(ncontact[i, iside], tncontact)
                np.copyto(fcontact[i, iside], tfcontact)
    # END: LOOP over frames

    if outtype == 'outb':  # Need to calculate ncontact & fcontact from un-normalized one
        ncontact, fcontact =\
            myvorn.update_contact_bilayer(
                framenum, ntype, ncontact, fcontact, tnmol_type)

    # Process to get statistics...
    if qa:
        ancontact = np.average(ncontact, axis=0)
        sncontact = np.std(ncontact, axis=0)

        afcontact = np.average(fcontact, axis=0)
        sfcontact = np.std(fcontact, axis=0)

        print('# Write average output')
        # write CONTACT outputs
        if outtype == 'outb':  # bilayer output
            obstype = 'ncomp'
            write_ave_std_bilayer(ntype, name_type, ancontact, sncontact,
                                  odir, obstype, suffix)

            obstype = 'fcomp'
            write_ave_std_bilayer(ntype, name_type, afcontact, sfcontact,
                                  odir, obstype, suffix)

        elif outtype == 'outl':  # leaflet output
            obstype = 'ncomp'
            write_ave_std_leaflet(nside, sside, ntype, name_type,
                                  ancontact, sncontact,
                                  odir, obstype, suffix)
            obstype = 'fcomp'
            write_ave_std_leaflet(nside, sside, ntype, name_type,
                                  afcontact, sfcontact,
                                  odir, obstype, suffix)

    # Time series output
    if qt:  # pass
        print('# Write time series output')
        # dt=1.0/float(framenum) # increment in time
        if outtype == 'outb':
            obstype = 'ncomp'
            write_time_series_bilayer(framenum, interval, time_step,
                                      ntype, name_type,
                                      ncontact, Snmol_type,
                                      odir, obstype, suffix)

            obstype = 'fcomp'
            write_time_series_bilayer(framenum, interval, time_step,
                                      ntype, name_type,
                                      fcontact, Snmol_type,
                                      odir, obstype, suffix)

        elif outtype == 'outl':
            obstype = 'ncomp'
            write_time_series_leaflet(framenum, interval, time_step,
                                      nside, sside, ntype, name_type,
                                      ncontact, SLnmol_type,
                                      odir, obstype, suffix)

            obstype = 'fcomp'
            write_time_series_leaflet(framenum, interval, time_step,
                                      nside, sside, ntype, name_type,
                                      fcontact, SLnmol_type,
                                      odir, obstype, suffix)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_voronoi_contact(**settings)


if __name__ == '__main__':
    main()
