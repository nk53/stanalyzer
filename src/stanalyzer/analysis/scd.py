#!/usr/bin/python
import argparse
import re
import typing as t

import numpy as np
import MDAnalysis as mda
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup

import stanalyzer.cli.stanalyzer as sta
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol
from . import scd_util as myscd
from .scd_util import StrList2D, StrList3D, IntList2D

ANALYSIS_NAME = 'scd'

if t.TYPE_CHECKING:
    import numpy.typing as npt

NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
NDInt64: t.TypeAlias = 'npt.NDArray[np.int64]'

LeafletAssignmentMethod: t.TypeAlias = t.Literal['mda', 'zpos']
OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']


# --- The following are hard set for membrane analysis
nside = 2     # up/dn
sside = ["up", "dn"]


class ProcessedArgs(t.NamedTuple):
    selection: list[str]
    ntype: int
    sel_type: list[str]
    name_type: list[str]
    ic_name: StrList2D
    nchain: list[int]
    qsplit: list[bool]


def process_args(sel: str, split_to_mol: str = '') -> ProcessedArgs:
    """
    ----------
    Process arguments
    ----------
    """

    # selections for individual lipid types
    selection: list[str] = re.split(';|,', f'{sel:s}')
    ntype = len(selection)                  # number of unique lipid types
    for i in range(0, ntype):
        selection[i] = selection[i].strip()

    # Process selection strings to extract information
    sel_type: list[str] = []        # selection type
    name_type: list[str] = []       # name of lipid type
    ic_name: StrList2D = []   # starting carbon names for ind. tail in ind. lipid type
    # number of carbon chains in ind. lipid types
    nchain: list[int] = []
    for i in range(0, ntype):
        tmps = selection[i].split()
        sel_type.append(tmps[0])   # segid/resname/moleculetype
        name_type.append(tmps[1])  # PROA/PRO*/DSPC/...
        ic_name.append([])         #

        # process remaining strings
        if sel_type[i] == 'segid':  # segid based selection
            if tmps[3] == 'resname':
                jst = 5  # segid XXX and resname YYY and (names)
            else:
                jst = 2  # segid XXX and (names)
        else:  # resname based selection
            jst = 2  # resname YYY ...
        for j in range(jst, len(tmps)):
            tmps_check = tmps[j].strip("(").strip(")")
            tstr = tmps_check.lower()
            if len(tstr) == 0:
                continue
            if tstr == "name" or tstr == "and" or tstr == "or" or tstr == "not":
                continue
            # if j < len(tmps) - 1:
            if j < len(tmps):
                ic_name[i].append(tmps_check)  # append starting carbon name
        # update nchain
        nchain.append(len(ic_name[i]))

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
    else:  # get args.split upto ntype
        for i in range(0, ntype):
            if split[i].lower() == "y":
                qsplit.append(True)
            else:
                qsplit.append(False)

    return ProcessedArgs(selection, ntype, sel_type, name_type, ic_name, nchain, qsplit)


def write_ave_std_leaflet(nside: int, sside: list[str],
                          ntype: int, name_type: list[str],
                          nchain: list[int], ncarbons: IntList2D, carbons: StrList3D,
                          array1: t.Any, array2: t.Any,
                          SLnmol_type: NDInt64,
                          odir: str, suffix: str) -> None:
    # sside : leafelt names
    # array1: average
    # array2: std
    for i in range(0, nside):
        side = sside[i]
        tnmol_type = SLnmol_type[i]

        # print for lipid types
        for j in range(0, ntype):
            tnamej = name_type[j].strip("*")

            # handel nmol_type = 0 case
            if tnmol_type[j] == 0:
                sout = f'no {tnamej} in {side} leaflet.\n'
                sta.write_to_outfile(
                    f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat',
                    sout)
                continue

            # Write output strings for nmol_type > 0
            ntail = nchain[j]
            for k in range(0, ntail):
                sout = f'# {side}: {tnamej}  chain{k} SCD\n'
                # loop over carbon
                nc = ncarbons[j][k]
                for m in range(0, nc):
                    sout += f' {carbons[j][k][m]:6s} {array1[i][j][k][m]:10.5f}' \
                        f' {array2[i][j][k][m]:10.5f}\n'
                # print(sout)
                sta.write_to_outfile(
                    f'{odir}/ave_{side}_{tnamej.lower()}_chain{k}_{suffix}.dat', sout)


def write_ave_std_bilayer(ntype: int, name_type: list[str],
                          nchain: list[int], ncarbons: IntList2D, carbons: StrList3D,
                          array1: t.Any, array2: t.Any,
                          Snmol_type: NDInt64,
                          odir: str, suffix: str) -> None:
    # array1: average
    # array2: std
    side = 'bilayer'
    # print for lipid types
    for j in range(0, ntype):
        tnamej = name_type[j].strip("*")

        # handel nmol_type = 0 case
        if Snmol_type[j] == 0:
            sout = f'no {tnamej} in the bilayer.\n'
            sta.write_to_outfile(
                f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat',
                sout)
            continue

        # Write output for nmol_type > 0
        ntail = nchain[j]
        for k in range(0, ntail):
            sout = f'# {side}: {tnamej} chain{k} SCD\n'
            # loop over carbon
            nc = ncarbons[j][k]
            for m in range(0, nc):
                sout += f' {carbons[j][k][m]:6s} {array1[j][k][m]:10.5f} {array2[j][k][m]:10.5f}\n'
            # print(sout)
            sta.write_to_outfile(
                f'{odir}/ave_{side}_{tnamej.lower()}_chain{k}_{suffix}.dat', sout)


def write_time_series_leaflet(framenum: int, interval: int, time_step: float,
                              nside: int, sside: list[str],
                              ntype: int, name_type: list[str],
                              nchain: list[int],
                              ncarbons: IntList2D, carbons: StrList3D,
                              array: t.Any,
                              SLnmol_type: NDInt64,
                              odir: str, suffix: str) -> None:
    # sside : leafelt names
    # array : time series
    for i in range(0, nside):
        tnmol_type = SLnmol_type[i]

        side = sside[i]
        for j in range(0, ntype):
            tnamej = name_type[j].strip("*")

            # handle nmol_type = 0 case
            if tnmol_type[j] == 0:
                sout = f'no {tnamej} in {side} leaflet.\n'
                sta.write_to_outfile(
                    f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat',
                    sout)
                continue

            # write output for nmol_type > 0
            ntail = nchain[j]
            for k in range(0, ntail):
                nc = ncarbons[j][k]

                # generate output string
                sout = f'# {side}: time series of {tnamej} SCD for chain {k}\n'
                sout += f'#{"frame":10s}'
                for n in range(0, nc):
                    sout += f' {carbons[j][k][n]:10}'
                sout += '\n'

                for m in range(0, framenum):
                    sout += f'{time_step*(interval*m+1):10.5f}'
                    # sout += f'{interval*m+1:10d}'
                    for n in range(0, nc):
                        sout += f' {array[i][j][k][n, m]:10.5f}'
                    sout += '\n'
                sta.write_to_outfile(
                    f'{odir}/time_{side}_{tnamej.lower()}_chain{k}_{suffix}.dat', sout)
                # print(sout)


def write_time_series_bilayer(framenum: int, interval: int, time_step: float,
                              ntype: int, name_type: list[str],
                              nchain: list[int],
                              ncarbons: IntList2D, carbons: StrList3D,
                              array: t.Any, Snmol_type: NDInt64,
                              odir: str, suffix: str) -> None:
    side = "bilayer"
    for j in range(0, ntype):
        tnamej = name_type[j].strip("*")

        # handle nmol_type = 0 case
        if Snmol_type[j] == 0:
            sout = f'no {tnamej} in the bilayer.\n'
            sta.write_to_outfile(
                f'{odir}/NA_{side}_{tnamej.lower()}_{suffix}.dat',
                sout)
            continue

        # write output for nmol_type > 0
        ntail = nchain[j]
        for k in range(0, ntail):
            nc = ncarbons[j][k]

            # generate output string
            sout = f'# {side}: time series of {tnamej} SCD for chain {k}\n'
            sout += 'f#{"time":10s}'
            for n in range(0, nc):
                sout += f' {carbons[j][k][n]:10s}'
            sout += '\n'

            for m in range(0, framenum):
                # sout += f'{interval*m}'
                for n in range(0, nc):
                    sout += f' {array[j][k][n, m]:10.5f}'
                sout += '\n'
            sta.write_to_outfile(
                f'{odir}/time_{side}_{tnamej.lower()}_chain{k}_{suffix}.dat', sout)
            # print(sout)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'interval', 'time_step')
    parser.add_argument('--sel', metavar='selection', default='',
                        help='Selection for individual molecule type with a '
                        'format, "segid/resname/moltype MOLECULENAME and name '
                        'CARBONNAMES". For the CHARMM format topology, a '
                        'selection for POPC can be "resname POPC and (name C22 '
                        'or name C32)".')
    parser.add_argument('--split', action='store',  # nargs='+',default='Y',
                        default='', help='Y/N. If N, the atom group for the selection '
                        'is considered as that for a single molecule. If Y,'
                        ' the atom group is further splitted to molecular '
                        'level based on segid/resname/moleculename. '
                        'Default is Y/y.')
    parser.add_argument('--sel-sys', metavar='selection', default='',
                        help='Selection for system atom groups in membranes'
                        ' for bilayer recentering.')
    # ## Commented out: Activate when non-planar bilayers are supported.
    # parser.add_argument('--lam', metavar='OPT', default='mda', choices=['mda', 'zpos'],
    #                     help='Leleat assignment method. mda: MDAnalysis.analysis.leaflet;'
    #                     ' zpos: Z-position. Default: mda')
    parser.add_argument('--suffix', type=str, default='0',
                        help='Suffix  to output file(s)')
    parser.add_argument('--otype', metavar='OPT', default='outl', choices=['outl', 'outb'],
                        help="Output type. outl: leaflets; outb: bilayer. Default: mda")
    parser.add_argument('--qa', action='store_true', default=False,
                        help='Write output average SCDs if provided')
    parser.add_argument('--qt', action='store_true', default=False,
                        help='Write output time series if provided')
    return parser


def run_scd(sel: str, split: str, qa: bool, qt: bool,
            sel_sys: str, psf: sta.FileRef, traj: sta.FileRefList,
            suffix: str, time_step: float | str,
            interval: int = 1,  # lam: LeafletAssignmentMethod='mda',
            otype: OutputFileType = 'outl') -> None:
    """
    ----------
    Calculate SCD order parameters.

    SCDs are calculated individual carbons in individual chains in individual molecules.
    Then, SCDs are avaraged for individual lipid types in individual leaflets.
    ----------
    """

    # process non-system arguments
    selection, ntype, sel_type, name_type, ic_name, nchain, qsplit =\
        process_args(sel, split)

    # method=lam
    method = 'zpos'  # use this option for planar bilaeyrs
    outtype = otype
    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    # print summary
    for i in range(0, ntype):
        print(f'sel: {selection[i]} '
              f' ; ic_name: {ic_name[i]} ; nchain = {nchain[i]}')
        print(f'#Split "{selection[i]}" into molecule level', qsplit[i])
    print(f'Bilayer is recentered at z = 0 using {sel_sys}:')
    if method == "zpos":
        print('Leaflets are assigned based on z-positions')
    elif method == "mda":
        print('Leaflets are assgined using a modified version of MDA LeafletFinder')
    print('Output type', outtype)
    print(f'SCD will be calculated every {interval} frames')
    print('Write average SCDs', qa)
    print('Write time series of SCD', qt)

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

    # Get numbers of molecules in individual lipid types,
    # Get total number of molecules,
    # Get lipid type index for individual molecules,
    # & Generate full atom groups for leaflet assignment
    #
    nmol_type, nmol, id_type, ag_full = \
        mymol.generate_full_mol_groups(u, ntype, sel_type, name_type, qsplit)

    # Generate
    #       carbon names for ind. chain in ind. lipid types
    #       hydrogen names boned to ind. carbons in ind. chain in ind. lipid types
    #       number of carbons in ind. chains in ind. lipid types
    #
    # carbons[type][chain][carbon]
    # hydrogens[type][chain][carbon][hydrogen]
    # ncarbons[type][chain]
    #
    carbons, hydrogens, ncarbons = \
        myscd.generate_carbons_hydrogens_ncarbons_type(u, ag_full, ntype, name_type,
                                                       nmol_type, ic_name, nchain)

    # Generate
    # 	atomgroups of carbons in ind. chains in ind. lipids
    # 	ag_c[lipid][chain]
    #   ag_h[lipid][chain][carbon]
    ag_c, ag_h = \
        myscd.generate_ag_carbon_hydrogen(u, ag_full, name_type, id_type,
                                          nchain, ncarbons, carbons, hydrogens)

    # raw SCD array for frames
    #  - raw data for individual carbon in indiviual tail in invidual lipids
    # raw_scd[lipid][chain][carbons] : list array upto lipid/chain level
    #                                : from carbons - numpy array
    #
    raw_scd = myscd.setup_raw_scd(nmol, id_type, nchain, ncarbons)

    # Set up output SCD
    scd, weight = myscd.setup_scd_weight(
        nside, ntype, nchain, ncarbons, framenum, outtype)
    ascd, astd = myscd.setup_ascd_astd(
        nside, ntype, nchain, ncarbons, outtype)  # average & std

    # Below, the membrane normal is assumed to be parallel to the z-axis
    memb_norm = np.array([0, 0, 1], dtype=float)

    # Analyze SCD
    # Loop over frames
    for i in range(0, framenum):
        print(f'# Analyze SCD: proceed {interval*i+1}/{interval*framenum}')
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

        # Calculate Raw SCD
        myscd.calculate_raw_scd(raw_scd, nmol, ag_c, ag_h, memb_norm)

        # Assign leaflet indices for individual lipids: Run this even qb = False
        if method == 'zpos':
            if i == 0:
                print('# leaflet assignemnt based on z-position')
            leaflets = myleaflet.assign_leaflet_zpos(u, ag_cent)
        elif method == 'mda':
            if i == 0:
                print(
                    '# leaflet assignemnt based on LeafletFinder with hydrid cutoff search')
            leaflets = myleaflet.assign_leaflet(u, ag_cent)
        id_side = mymol.assign_leaflet_index_to_full_ag(ag_full, leaflets)

        # Now process to get scd for individual lipid types for the frame, i
        scd, weight = myscd.calculate_scd(
            i, ntype, raw_scd, scd, weight, id_side, id_type, nchain, ncarbons, outtype)

    # weight[type][(leaflet,) frame]
    #     => Total sum: Snmol_type (bilayer)/SLnmol_type(leaflets)
    if outtype == 'outb':
        Snmol_type = np.zeros([ntype], dtype=int)
        for itype in range(0, ntype):
            Snmol_type[itype] = np.sum(weight[itype])
    elif outtype == 'outl':
        SLnmol_type = np.zeros([nside, ntype], dtype=int)
        for iside in range(0, nside):
            for itype in range(0, ntype):
                SLnmol_type[iside][itype] = np.sum(weight[iside][itype])

    # Calculate weighted ave and std of SCD
    ascd, astd = myscd.calculate_ave_and_std(
        ascd, astd, scd, weight, framenum, outtype)

    # Write output (for each lipid type)
    if outtype == 'outb':
        if qa:
            print('# Write average SCD for bilayer')
            write_ave_std_bilayer(ntype, name_type, nchain, ncarbons, carbons,
                                  ascd, astd, Snmol_type, odir, suffix)
        if qt:
            print('# Write time series output')
            write_time_series_bilayer(framenum, interval, time_step,
                                      ntype, name_type,
                                      nchain, ncarbons, carbons,
                                      scd, Snmol_type, odir, suffix)

    elif outtype == 'outl':
        if qa:
            print('# Write average SCD for leaflets')
            write_ave_std_leaflet(nside, sside, ntype, name_type,
                                  nchain, ncarbons, carbons,
                                  ascd, astd, SLnmol_type, odir, suffix)
        if qt:
            print('# Write time series output')
            write_time_series_leaflet(framenum, interval, time_step,
                                      nside, sside, ntype, name_type,
                                      nchain, ncarbons, carbons,
                                      scd, SLnmol_type, odir, suffix)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_scd(**settings)


if __name__ == '__main__':
    main()
