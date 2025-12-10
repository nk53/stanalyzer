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

if t.TYPE_CHECKING:
    import numpy.typing as npt

ANALYSIS_NAME = 'thickness'

# LeafletAssignmentMethod: t.TypeAlias = t.Literal['mda', 'zpos']
NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'

# --- The following are hard set for membrane analysis
nside = 2  # up/dn
sside = ["up", "dn"]


class ProcessedArgSys(t.NamedTuple):
    selection: list[str]
    ntype: int
    qsplit: list[bool]
    sel_type: list[str]
    name_type: list[str]


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


def write_output(out: sta.FileRef, framenum: int, interval: int, time_step: float,
                 zup: NDFloat64, zdn: NDFloat64, thick: NDFloat64,
                 azup: float, szup: float, azdn: float, szdn: float,
                 athick: float, sthick: float) -> None:

    sout = f'# Ave. zup = {azup:10.5f} (STD = {szup:10.5f})\n'
    sout += f'# Ave. zdn = {azdn:10.5f} (STD = {szdn:10.5f})\n'
    sout += f'# Ave. thick = {athick:10.5f}  STD = {sthick:10.5f}\n'
    sout += '#\n'
    sout += f'# {"time":8s} {"z_up":10s} {"z_dn":10s} {"thickness":10s}\n'
    for i in range(0, framenum):
        sout += f'{time_step*(interval*i+1):10.5f}'
        # sout += f'{interval*i+1:10d}'
        sout += f' {zup[i]:10.5f} {zdn[i]:10.5f} {thick[i]:10.5f}\n'

    with sta.resolve_file(out, 'w') as outfile:
        print(sout, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out', 'interval', 'time_step')
    parser.add_argument('--sel', metavar='selection', default='',
                        help="Membrane headgroup atom selection")
    parser.add_argument('--sel-sys', metavar='selection', default='',
                        help='Selection for system atom groups for leaflet assignment'
                        'and centering.')
    # ## Commented out: Activate when non-planar bilayers are supported.
    # parser.add_argument('--lam', metavar='OPT', default='mda', choices=['mda', 'zpos'],
    #                     help='Leleat assignment method. mda: MDAnalysis.analysis.leaflet; '
    #                     'zpos: Z-position. Default: mda')
    parser.add_argument('--qcent', default=False, action='store_true',
                        help='If set True, bilayer midplane will be centered '
                        'using leaflet z-density profiles.')

    return parser


def run_thickness(sel: str, sel_sys: str, qcent: bool,
                  psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                  time_step: float | str,
                  # lam: LeafletAssignmentMethod='mda',
                  interval: int = 1,
                  ) -> None:
    """
    ----------
    Calculate Thickness of selected atom groups

    Leaflet thicknesses are calculated as mean z-positions of leaflet atom groups
    Bilayer thickness is the difference of two.
    ----------
    """

    # process sys arguments - To generate full ag groups for individual molecules
    selection_sys, ntype_sys, qsplit_sys, sel_type_sys, name_type_sys \
        = process_arg_sys(sel_sys)

    # method=lam
    method = 'zpos'  # use this option for plana bilayers
    if isinstance(time_step, str):
        time_step = float(time_step.split()[0])

    # READ topology and trajectory
    u = mda.Universe(psf, traj)  # MDA universe
    # number of frames to be analyzed
    framenum = int(u.trajectory.n_frames/interval)

    ag_sel = u.select_atoms(f'{sel:s}')  # selection for thickness analysis

    # centering
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

    if qcent:
        print('### Generation of full atom groups for membrane molecules: START')
        # Get numbers of molecules in individual lipid types,
        # Get total number of molecules,
        # Get lipid type index for individual molecules,
        # & Generate full atom groups for leaflet COM drift correction
        #
        nmol_type_sys, nmol_sys, id_type_sys, ag_full_sys = \
            mymol.generate_full_mol_groups(
                u, ntype_sys, sel_type_sys, name_type_sys, qsplit_sys)

        print('### Generation of full atom groups for membrane molecules: DONE')

    # output zup/zdn,thickness
    zup = np.zeros([framenum], dtype=float)
    zdn = np.zeros([framenum], dtype=float)
    thick = np.zeros([framenum], dtype=float)

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
        # xtla, xtlb, xtlc = ts.dimensions[:3]
        # box = np.array([xtla, xtlb, xtlc], dtype=float)

        # select up/dn leaflet of sel
        if method == 'zpos':
            ag_leaflet = myleaflet.assign_leaflet_zpos(u, ag_sel)
        elif method == 'mda':
            ag_leaflet = myleaflet.assign_leaflet(u, ag_sel)

        # thicknesses
        _, _, tzup = ag_leaflet[0].center_of_geometry()
        _, _, tzdn = ag_leaflet[1].center_of_geometry()
        zup[i] = tzup
        zdn[i] = tzdn
        thick[i] = tzup - tzdn

    # ave & std of zup/zdn/thickness
    azup = np.average(zup)
    szup = np.std(zup)
    azdn = np.average(zdn)
    szdn = np.std(zdn)
    athick = np.average(thick)
    sthick = np.std(thick)

    # write output
    write_output(out, framenum, interval, time_step,
                 zup, zdn, thick,
                 azup, szup, azdn, szdn, athick, sthick)


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginning of this function
    run_thickness(**settings)


if __name__ == '__main__':
    main()
