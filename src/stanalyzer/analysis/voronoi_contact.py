import argparse
import os
import re

import MDAnalysis as mda
import MDAnalysis.transformations as transformations
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from scipy.spatial import Voronoi

import stanalyzer.cli.stanalyzer as sta
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol
from . import voronoi_analysis as myvorn

ANALYSIS_NAME = 'voronoi_apl'


# --- The following are hard set for membrane analysis
dim = 2              # 2 dimension
nimage = int(3**dim)  # number of total primary+images for 2D

nside = 2     # up/dn
sside = ["up", "dn"]


def process_args(sel, split):
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

    split = re.split(';|,', f'{split:s}')
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

    return selection, ntype, qsplit


def write_ave_std_leaflet(nside, ntype, sside, name_type, array1, array2, odir, otype):
    # sside:  leaflet name
    # array1: average
    # array2: std
    for i in range(0, nside):
        side = sside[i]
        sout = f'# {side}: {otype} of contacting lipid types & std\n'
        for j in range(0, ntype):
            sout += f' {name_type[j]:21s}'
        sout += '\n'
        for j in range(0, ntype):
            for k in range(0, ntype):
                sout += f' {array1[i, j, k]:10.5f} {array2[i, j, k]:10.5f}'
            sout += '\n'

        fout = f'{odir}/{side}_{otype}.plo'
        f = open(fout, 'w')
        f.write(sout)
        f.close()
        print(sout)


def write_ave_std_bilayer(ntype, name_type, array1, array2, odir, otype):
    # array1: average
    # array2: std
    side = "bilayer"
    sout = f'# {side}: {otype} of contacting lipid types & std\n'
    for j in range(0, ntype):
        sout += f' {name_type[j]:21s}'
    sout += '\n'
    for j in range(0, ntype):
        for k in range(0, ntype):
            sout += f' {array1[j, k]:10.5f} {array2[j, k]:10.5f}'
        sout += '\n'

    fout = f'{odir}/{side}_{otype}.plo'
    f = open(fout, 'w')
    f.write(sout)
    f.close()
    print(sout)


def write_time_series_leaflet(framenum, interval, nside, ntype,
                              sside, name_type, array, odir, otype):
    # array: time series
    for j in range(0, nside):
        side = sside[j]
        for k in range(0, ntype):
            sout = f'# {side}: {otype} of contacting lipid types for {name_type[k]}\n'
            sout += '#     frame'
            for m in range(0, ntype):
                sout += f' {name_type[m]:10s}'
            sout += '\n'
            for i in range(0, framenum):
                # ct = (dt*(i+1)+(cnt-1))
                # sout += f' {ct}'
                sout += f' {interval*i:10d}'
                for m in range(0, ntype):
                    sout += f' {array[i, j, k, m]:10.5f}'
                sout += '\n'

            fout = f'{odir}/time_{side}_{name_type[k]}_{otype}.plo'
            f = open(fout, 'w')
            f.write(sout)
            f.close()
            # print(sout)


def write_time_series_bilayer(framenum, interval, ntype, name_type, array, odir, otype):
    # array: time series
    side = "bilayer"
    for j in range(0, ntype):
        sout = f'# {side}: {otype} of contacting lipid types for {name_type[j]}\n'
        sout += '#     frame'
        for k in range(0, ntype):
            sout += f' {name_type[k]:10s}'
        sout += '\n'
        for i in range(0, framenum):
            # ct = (dt*(i+1)+(cnt-1))
            # sout += f' {ct}'
            sout += f' {interval*i:10d}'
            for k in range(0, ntype):
                sout += f' {array[i, j, k]:10.5f}'
            sout += '\n'

        fout = f'{odir}/time_{side}_{name_type[j]}_{otype}.plo'
        f = open(fout, 'w')
        f.write(sout)
        f.close()
        # print(sout)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval')
    parser.add_argument('--center', default=False, action='store_true',
                        help='Perform membrane centering.')

    parser.add_argument('--sel', metavar='selection',
                        help='Selection for individual molecule type with a format, '
                        'segid/resname/moltype MOLECULENAME and name ATOMNAMES. For the CHARMM '
                        'format topology, an selection for POPC can be "resname POPC and (name C2 '
                        'or name C21 or name C31)".')
    parser.add_argument('--split', action='store',
                        help='Y/N. If N, the atom group for the selection is considered as that '
                        'for a single molecule. If Y, the atom group is further splitted to '
                        'molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--sel-sys', metavar='selection',
                        help='Selection for system atom groups in membranes for bilayer '
                        'recentering.')
    parser.add_argument('--qz', action='store_true', default=False,
                        help='Z-position based leaflet assigment: Default is false. Maybe useful '
                        'when selections are minor components of the bilayer')
    parser.add_argument('--qb', action='store_true', default=False,
                        help='Do bilayer analysis if provided: only for sym. bilayers')
    parser.add_argument('--qt', action='store_true', default=False,
                        help='Write output time series if provided')
    parser.add_argument('--qa', action='store_true', default=False,
                        help='Write output averages if provided')

    return parser


def run_voronoi_contact(sel, split, qz, qb, qt, qa, sel_sys, psf: sta.FileRef, traj:
                        sta.FileRefList, interval: int = 1, center: bool = False) -> None:
    """
    ----------
    Analyze contacts between inidividual molecul types.

    Contacts are calculated between individual molecules in individual leaflets.
    Then, avaraged for contact between individual lipid types.
    ----------
    """

    # process non-system arguments
    selection, ntype, qsplit = process_args(sel, split)

    # print summary of arguments
    for i in range(0, ntype):
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    print('Write results for bilayer', qb)
    print('Write averge contact numbers/fractions', qa)
    print('Write time sereis output', qt)
    print('Bilayer is recentered at z = 0 using {sel_sys}:', center)
    print('Analysis will be done every {interval}frames')

    # make output dir
    odir = "./voronoi/contact"

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    # if center: bilayer recentering - should be done before any assignments
    # - center in the box (an atom)
    # - center in the box (atom group for system)
    # - unwrap to get connectd molecules
    # Taken from
    if center:
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
    if qb:
        ncontact = np.zeros([framenum, ntype, ntype], dtype=float)
        fcontact = np.zeros([framenum, ntype, ntype], dtype=float)
        # normalization factor
        tnmol_type = np.zeros([framenum, ntype], dtype=float)
    else:  # if (qb is False):
        ncontact = np.zeros([framenum, nside, ntype, ntype], dtype=float)
        fcontact = np.zeros([framenum, nside, ntype, ntype], dtype=float)

    # START: LOOP over frames
    for i in range(0, framenum):
        #  ct=(cnt-1)+dt*(i+1) # in ns
        print(f'# processing {interval*i+1}/{interval*framenum}')
        ts = u.trajectory[interval*i]

        # do frame-wise bilayer recentering - remaining translation
        if center:
            Lag_ref = myleaflet.assign_leaflet_zpos(u, ag_cent)
            zref = np.zeros([2], dtype=float)
            for i in range(0, nside):
                zref[i] = np.mean(Lag_ref[i].positions[:, 2])
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
                u, nside, ntype, selection, qsplit, sside, qz)

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
                sinfo += f' {name_type[itype]}'
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
            if qb:
                ncontact[i] += (tncontact.T * nmol_type).T
                tnmol_type[i] += nmol_type
            else:  # if (qb is False):
                np.copyto(ncontact[i, iside], tncontact)
                np.copyto(fcontact[i, iside], tfcontact)
    # END: LOOP over frames

    if qb:  # Need to calculate ncontact & fcontact from un-normalized one
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
        if qb:
            otype = "numb"
            write_ave_std_bilayer(
                ntype, name_type, ancontact, sncontact, odir, otype)

            otype = "frac"
            write_ave_std_bilayer(
                ntype, name_type, afcontact, sfcontact, odir, otype)
        else:  # if (qb is False):
            otype = "numb"
            write_ave_std_leaflet(nside, ntype, sside,
                                  name_type, ancontact, sncontact, odir, otype)

            otype = "frac"
            write_ave_std_leaflet(nside, ntype, sside,
                                  name_type, afcontact, sfcontact, odir, otype)

    # Time series output
    if qt:  # pass
        print('# Write time series output')
        # dt=1.0/float(framenum) # increment in time
        if qb:
            otype = "numb"
            write_time_series_bilayer(
                framenum, interval, ntype, name_type, ncontact, odir, otype)

            otype = "frac"
            write_time_series_bilayer(
                framenum, interval, ntype, name_type, fcontact, odir, otype)
        else:  # if (qb is False):
            otype = "numb"
            write_time_series_leaflet(
                framenum, interval, nside, ntype, sside, name_type, ncontact, odir, otype)

            otype = "frac"
            write_time_series_leaflet(
                framenum, interval, nside, ntype, sside, name_type, fcontact, odir, otype)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_voronoi_contact(**settings)


if __name__ == '__main__':
    main()
