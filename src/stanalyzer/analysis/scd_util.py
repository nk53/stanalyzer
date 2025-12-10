#!/usr/bin/python
import sys
import math
import re
import numpy as np

import typing as t
if t.TYPE_CHECKING:
    import numpy.typing as npt
    from MDAnalysis import AtomGroup, Universe

T = t.TypeVar('T')
Tup2: t.TypeAlias = tuple[T, T]
List2D: t.TypeAlias = list[list[T]]
List3D: t.TypeAlias = list[list[list[T]]]
List4D: t.TypeAlias = list[list[list[list[T]]]]
StrList2D: t.TypeAlias = List2D[str]
StrList3D: t.TypeAlias = List3D[str]
StrList4D: t.TypeAlias = List4D[str]
IntList2D: t.TypeAlias = List2D[int]
IntList3D: t.TypeAlias = List3D[int]
AGroup2D: t.TypeAlias = list[list['AtomGroup']]
AGroup3D: t.TypeAlias = list[list[list['AtomGroup']]]
NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
SmallWeightTup: t.TypeAlias = tuple[List2D[NDFloat64], list[NDFloat64]]
BigWeightTup: t.TypeAlias = tuple[List3D[NDFloat64], List2D[NDFloat64]]
WeightTup: t.TypeAlias = SmallWeightTup | BigWeightTup
List2DTup: t.TypeAlias = Tup2[List2D[NDFloat64]]
List3DTup: t.TypeAlias = Tup2[List3D[NDFloat64]]

OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']


class CHInner(t.NamedTuple):
    carbons: StrList2D
    hydrogens: StrList3D  # StrList2D
    ncarbons: list[int]


class CHOuter(t.NamedTuple):
    carbons: StrList3D
    hydrogens: StrList4D  # StrList3D
    ncarbons: IntList2D


# nside=2     # up/dn
def find_carbons_hydrogens(u: 'Universe', atomgroup: 'AtomGroup',
                           name_t: str, c_st: list[str], ntail: int) -> CHInner:
    """
    ----------
    Find carbons and associated bonded hydrogens in tails of a given molecule
    ----------

    input
          u        : MDA Universe
          atomgroup: full atom group for molecule
          name_t   : molecule name
          c_st     : starting carbon name for individual chains in the lipid
          ntail    : number of tails in the lipid type

    output
          carbons  : a list array of carbon names in invidual chains
          hydrogens: a list array of hydrogen names bonded to carbon in individual chains
          ncarbons : a list array of number of carbon atoms in individual chains
    """

    # list arrays of
    # carbons and hydrogens for the given lipid type
    carbons: StrList2D = []     # [nchain][0-ncarbons]
    hydrogens: StrList3D = []   # [nchain][0-ncarbons][...]
    ncarbons: list[int] = []

    # expand arrays along chains
    for i in range(0, ntail):
        carbons.append([])
        hydrogens.append([])
        ncarbons.append(0)

    # Search for carbons
    for i in range(0, ntail):
        c_curr = c_st[i]  # starting carbon name in the chain

        # start from ic_name of each chain
        # until there exist no more connected carbon atoms to the c_curr
        tag_c = atomgroup.intersection(u.select_atoms(f'name {c_curr}'))
        tag_h = atomgroup.intersection(
            u.select_atoms(f'name H* and bonded name {c_curr}'))

        if len(tag_c) == 0:
            print(f'# Starting Carbon name {c_curr} does not exist!')
            sys.exit(1)
        if len(tag_h) == 0:
            # sys.exit(1)
            print(f'# There exist no hydrogens bonded to {c_curr}')

        carbons[i].append(c_curr)        # first carbon atom in the chain
        hydrogens[i].append(tag_h.names)  # hydrogens bonded to carbon, c_curr
        ncarbons[i] += 1                 # update number of carbons

        # determine how to add more carbons
        # This is not unique:
        # e.g., C22,C23,... or C4S,C5S,...
        #
        # qsuffix = True
        if re.search(r'[0-9]+', c_curr[-1]):
            qsuffix = True
        else:
            qsuffix = False

        if qsuffix:
            icnt_c = int(c_curr[-1])     # starting carbon number
            cprefix = c_curr[0:-1]       # prefix for carbon name
            csuffix = ""                 # suffix
        else:
            icnt_c = int(c_curr[1:-1])  # starting counter for carbon atoms
            cprefix = c_curr[0]
            csuffix = c_curr[-1]
        # print(c_st,icnt_c) ; sys.exit(0);
        while True:
            c_prev = c_curr
            c_curr = f'{cprefix}{icnt_c+1}{csuffix}'
            tag_c = atomgroup.intersection(
                u.select_atoms(f'name {c_curr} and bonded name {c_prev}'))
            tag_h = atomgroup.intersection(
                u.select_atoms(f'name H* and bonded name {c_curr}'))
            # print(tag_c.names,len(tag))
            if len(tag_c) > 0:
                carbons[i].append(c_curr)        # append new carbon
                # append new hydrogens bonded to carbon, c_curr
                hydrogens[i].append(tag_h.names)
                ncarbons[i] += 1                 # update the number of carbons
                icnt_c += 1
            else:
                tname = name_t.strip("*")
                print(f'# {tname} chain {i}: all carbon atoms assigned')
                break

        # Print out assigned carbon/bonded hydrogens
        for j in range(0, ncarbons[i]):
            tname = name_t.strip("*")
            print(f'# {tname}: chain {i} ica = {j}',
                  carbons[i][j], hydrogens[i][j])

    return CHInner(carbons, hydrogens, ncarbons)


def generate_carbons_hydrogens_ncarbons_type(
        u: 'Universe', ag_full: list['AtomGroup'], ntype: int,
        name_type: list[str], nmol_type: list[int], ic_name: StrList2D,
        nchain: list[int]) -> CHOuter:
    """
    ----------
    Generate carbons,hydroges,ncarbons for ind. chain in ind. lipid types
    Find carbons in individual chains in individual lipid types
    These will be used to generate atom groups for carbons and bonded hydrogens
    ----------

    input
          u      : MDA Universe
          ag_full: atom group of a molecule with all atoms for all moleules
          ntype  : numbrer of molecule types
          name_type: names of individual molecule types
          nmol_type: numbers of molecules of individual molecule types
          ic_name  : names of starting carbons in individual chains in ind. mol. types
          nchain   : numbers of chains in individual molecule types

    output
          carbons  : carbon names in ind. chains in ind. molecule types
          hydrogens: hydrogen names boneded to carbons
          ncarbons : number of carbons in ind. chains in ind. molecule types
    """

    carbons: StrList3D = []    # [ntype][nchain][0-ncarbons]
    hydrogens: StrList4D = []  # [ntype][nchain][0-ncarbons][...]
    ncarbons: IntList2D = []

    imol = 0
    for i in range(0, ntype):
        carbons.append([])
        hydrogens.append([])
        ncarbons.append([])

        # expand arrays along chains
        ntail = nchain[i]
        carbons[i], hydrogens[i], ncarbons[i] = find_carbons_hydrogens(
            u, ag_full[imol], name_type[i], ic_name[i], ntail)

        imol += nmol_type[i]

    return CHOuter(carbons, hydrogens, ncarbons)


def generate_ag_carbon_hydrogen(
        u: 'Universe', atomgroups: list['AtomGroup'], name_type:
        list[str], id_type: list[int], nchain: list[int],
        ncarbons: IntList2D, carbons: StrList3D,
        hydrogens: StrList4D) -> tuple['AtomGroup', 'AtomGroup']:
    """
    ----------
    Generate atom groups of carbon/associated bonded hydrogens for individual lipids
    ----------

    input
          u         : MDA Universe
          atomgroups: atomgroups of individual molecules
          nchain    : number of chains in individual molecule types
          ncarbons  : number of carbons in ind. chains in unique mol. types
          carbons   : carbon names in ind. chain in unique molecule types
          hydrogens : hydrogen names bonded to ind. carbons in ind. chains

    output
          ag_c      : atom groups of carbons in ind. chain in ind. lipids
          ag_h      : atom groups of Hs bonded to ind. Cs in ind. tail in ind. lipids
    """

    nmol = len(atomgroups)

    # atom group list arrays
    ag_c: AGroup2D = []
    ag_h: AGroup3D = []

    # Loop over lidids
    for i in range(0, nmol):
        ag_c.append([])
        ag_h.append([])

        itype = id_type[i]                    # lipid type index
        name_t = name_type[itype].strip("*")  # molecule type name

        ntail = nchain[itype]  # number of chains in the lipid, i
        if i % int(nmol/10) == 0:
            print(f'# GENRATION OF ATOM GROUP for {i}/{nmol}: {name_t}')
        # loop over chains
        for j in range(0, ntail):
            ag_c[i].append([])
            ag_h[i].append([])

            ag_c[i][j] = u.atoms[[]]  # an empty atomgroup

            c_names = carbons[itype][j]  # carbon names
            nc = ncarbons[itype][j]     # number of carbons
            # loop over carbon atoms
            for k in range(0, nc):
                tc = c_names[k]
                sel_c = f'name {tc}'
                tag_c = atomgroups[i].intersection(u.select_atoms(sel_c))

                h_names = hydrogens[itype][j][k]
                nh = len(h_names)
                sel_h = ""
                for m in range(0, nh):
                    if m == 0:
                        sel_h += f'name {h_names[m]}'
                    else:
                        sel_h += f' or name {h_names[m]}'
                tag_h = atomgroups[i].intersection(u.select_atoms(sel_h))

                ag_c[i][j] += tag_c      # adding carbons in the chain
                # appending bonded hydrogens to the carbon
                ag_h[i][j].append(tag_h)

    return ag_c, ag_h


def setup_raw_scd(nmol: int, id_type: list[int], nchain: list[int],
                  ncarbons: IntList2D) -> List2D[NDFloat64]:
    """
    ----------
    Set up raw SCD array
    ----------
    NOTE:
    raw SCD array
     - raw data for individual carbon in indiviual tails in invidual lipids
    raw_scd[lipid][chain][carbons] : list array upto lipid/chain level
                                   : from carbons - numpy array

    input
          nmol    : number of lipids
          id_type : list array of the lipid type of individual lipids
          nchain  : list array of the number of chains array for individual lipids
          ncarbons: list array of the number of carbons in individual chains

    output
          raw_scd : list array of raw scd
    """

    raw_scd: List2D[NDFloat64] = []
    for i in range(0, nmol):
        # lipid type and associated number of tails
        itype = id_type[i]
        ntail = nchain[itype]

        raw_scd.append([])
        for j in range(0, ntail):
            nc = ncarbons[itype][j]

            # now make numpy array ## raw_scd[i][j][0-ncarbons]
            raw_scd[i].append(np.zeros([nc], dtype=float))
    return raw_scd


def calculate_raw_scd(raw_scd: List2D[NDFloat64], nmol: int,
                      ag_c: List2D['AtomGroup'], ag_h: List3D['AtomGroup'],
                      memb_norm: NDFloat64) -> None:
    """
    ----------
    Calculate raw SCD for ind. carbons in ind. chains in ind. molecules
    ----------

    input
          nmol     : number of lipids
          ag_c     : Carbon atomgroups
          ag_h     : Hydrogen atomgroup associated to carbons
          memb_norm: membrane normal vector

    input/output
          raw_scd : SCD array for ind. Cs in ind. tails in ind. molecules
    """

    # Loop over lipids
    for j in range(0, nmol):
        # for j in range(0,1):
        # Loop over tails
        ntail = len(ag_c[j])  # number of tails
        for k in range(0, ntail):

            pos_c = ag_c[j][k].positions  # carbon position

            # Loop over carbons
            nc = len(ag_c[j][k])
            for m in range(0, nc):
                # CD vector
                vec = ag_h[j][k][m].positions  # hydrogen positions
                ref = pos_c[m]                # carbon position
                vec -= ref                    # CD vector

                # SCD = 1/2 * ( 3 * cos^2(thet) - 1)
                # angle between the membrane normal and CD vector
                cos_thet2 = np.inner(vec, memb_norm) / \
                    np.linalg.norm(vec, axis=1)
                cos_thet2 = np.square(cos_thet2)
                tscd = 0.5 * (3.0*np.mean(cos_thet2) - 1.0)

                # put SCD
                # raw_scd[j][k][m][i]=tscd
                raw_scd[j][k][m] = tscd


@t.overload
def setup_scd_weight(nside: int, ntype: int, nchain: list[int],
                     ncarbons: IntList2D, framenum: int,
                     otype: t.Literal['outb']) -> SmallWeightTup: ...


@t.overload
def setup_scd_weight(nside: int, ntype: int, nchain: list[int],
                     ncarbons: IntList2D, framenum: int,
                     otype: t.Literal['outl']) -> BigWeightTup: ...


def setup_scd_weight(nside: int, ntype: int, nchain: list[int],
                     ncarbons: IntList2D, framenum: int,
                     otype: OutputFileType) -> WeightTup:
    """
    ----------
    Set up SCD time series arrays and associated un-normalized weights
    ----------

    input
          nside   : number of leaflet
          ntype   : number of unique lipid types
          nchain  : number of chains in unique lipid types
          ncarbons: number of carbons in ind. chain in unique lipid types
          otype : output type (outl=Leaflets)/(outb=Bilayer)

    output
          scd     : SCD of ind. carbon in ind. chain in unique lipid types (bilayer/leaflet)
          weight  : number of counts for unique lipid types (bilayer/leaflet)
    """
    def do_small() -> SmallWeightTup:
        scd: List2D[NDFloat64] = []
        weight: list[NDFloat64] = []
        for i in range(0, ntype):
            scd.append([])
            weight.append(np.zeros([framenum], dtype=float))

            ntail = nchain[i]
            for j in range(0, ntail):
                nc = ncarbons[i][j]
                scd[i].append(np.zeros([nc, framenum], dtype=float))
        return scd, weight

    def do_big() -> BigWeightTup:
        scd: List3D[NDFloat64] = []
        weight: List2D[NDFloat64] = []
        for i in range(0, nside):
            scd.append([])
            weight.append([])
            for j in range(0, ntype):
                scd[i].append([])
                weight[i].append(np.zeros([framenum], dtype=float))

                ntail = nchain[j]
                for k in range(0, ntail):
                    nc = ncarbons[j][k]
                    scd[i][j].append(np.zeros([nc, framenum], dtype=float))
        return scd, weight

    if otype == 'outb':
        return do_small()
    return do_big()


@t.overload
def setup_ascd_astd(nside: int, ntype: int, nchain: list[int],
                    ncarbons: IntList2D, otype: t.Literal['outb']) -> List2DTup: ...


@t.overload
def setup_ascd_astd(nside: int, ntype: int, nchain: list[int],
                    ncarbons: IntList2D, otype: t.Literal['outl']) -> List3DTup: ...


def setup_ascd_astd(nside: int, ntype: int, nchain: list[int],
                    ncarbons: IntList2D, otype: OutputFileType) -> List2DTup | List3DTup:
    """
    ----------
    Setup average and std of SCDs
    ----------

    input
          nside   : number of leaflets
          ntype   : number of unique lipid types
          nchain  : number of chains in unique lipid types
          ncarbons: number of carbons in ind. chains in ind. lipid types
          otype : output type (outl=Leaflets)/(outb=Bilayer)

    output
          ascd    : average SCD
          astd    : std of SCD
    """
    def do_small() -> List2DTup:
        ascd: List2D[NDFloat64] = []  # weighted mean
        astd: List2D[NDFloat64] = []  # weighted std
        for i in range(0, ntype):
            ascd.append([])
            astd.append([])

            ntail = nchain[i]
            for j in range(0, ntail):
                nc = ncarbons[i][j]
                ascd[i].append(np.zeros([nc], dtype=float))
                astd[i].append(np.zeros([nc], dtype=float))
        return ascd, astd

    def do_big() -> List3DTup:
        ascd: List3D[NDFloat64] = []  # weighted mean
        astd: List3D[NDFloat64] = []  # weighted std
        for i in range(0, nside):
            ascd.append([])
            astd.append([])
            for j in range(0, ntype):
                ascd[i].append([])
                astd[i].append([])

                ntail = nchain[j]
                for k in range(0, ntail):
                    nc = ncarbons[j][k]
                    ascd[i][j].append(np.zeros([nc], dtype=float))
                    astd[i][j].append(np.zeros([nc], dtype=float))
        return ascd, astd

    if otype == 'outb':
        return do_small()
    return do_big()


def calculate_scd(iframe: int, ntype: int, raw_scd: List2D[NDFloat64],
                  scd: t.Any, weight: t.Any,
                  id_side: list[int], id_type: list[int], nchain: list[int],
                  ncarbons: IntList2D, otype: OutputFileType) -> WeightTup:
    """
    ----------
    Calculate SCD from raw SCD
    ----------
    input
          iframe  : frame index
          ntype   : number of unique lipid types
          raw_scd : frame raw SCD of individual lipids
          scd     : SCD for individual lipid types
          weight  : counts of individual lipid types
          id_side : leaflet indices of individual lipids
          id_type : lipid type inices of individual lipids
          nchain  : number of chains in individual lipid types
          ncarbons: number of carbons in individual chain for individual lipid types
          otype : output type (outl=Leaflets)/(outb=Bilayer)

    output
          scd,weight
    """

    nmol = len(raw_scd)
    # ntype: global parameter
    # nsdie: global parameter

    def do_small(scd: List2D[NDFloat64], weight: list[NDFloat64]) -> SmallWeightTup:
        for i in range(0, nmol):
            ii = id_type[i]
            ntail = nchain[ii]
            weight[ii][iframe] += 1.0

            for j in range(0, ntail):
                nc = ncarbons[ii][j]
                for k in range(0, nc):
                    scd[ii][j][k][iframe] += raw_scd[i][j][k]

        # normalize
        for i in range(0, ntype):
            ntail = nchain[i]
            tnorm = weight[i][iframe]

            if tnorm > 0.0:
                for j in range(0, ntail):
                    scd[i][j][:, iframe] /= tnorm
        return scd, weight

    def do_big(scd: List3D[NDFloat64], weight: List2D[NDFloat64]) -> BigWeightTup:
        nside = len(scd)
        for i in range(0, nmol):
            ii = id_type[i]
            iside = id_side[i]
            ntail = nchain[ii]
            weight[iside][ii][iframe] += 1.0
            for j in range(0, ntail):
                nc = ncarbons[ii][j]
                for k in range(0, nc):
                    scd[iside][ii][j][k][iframe] += raw_scd[i][j][k]

        # normalize to get ave
        for i in range(0, nside):
            for j in range(0, ntype):
                ntail = nchain[j]
                tnorm = weight[i][j][iframe]

                if tnorm > 0.0:
                    for k in range(0, ntail):
                        scd[i][j][k][:, iframe] /= tnorm
        return scd, weight

    if otype == 'outb':
        return do_small(
            t.cast(List2D[NDFloat64], scd),
            t.cast(list[NDFloat64], weight))

    return do_big(
        t.cast(List3D[NDFloat64], scd),
        t.cast(List2D[NDFloat64], weight))


def calculate_ave_and_std(ascd: t.Any, astd: t.Any, scd: t.Any, weight: t.Any,
                          framenum: int, otype: OutputFileType) -> tuple[t.Any, t.Any]:
    """
    ----------
    Calculate average and std of SCDs for individaul lipid types
    ----------
    input
          scd   : time series of SCD for ind. carbons in ind. chains in unique lipid types
          weight: un-normalized wight for each time frame
          otype : output type (outl=Leaflets)/(outb=Bilayer)

    input/output
          ascd  : weighted ave SCD for ind. carbons in ind. chains in unique lipid types
          astd  : weighted std of SCD
    """
    # nside/ntype : global parameters
    if otype == 'outb':
        weight = t.cast(list[NDFloat64], weight)
        # normalize weight
        ntype = len(weight)
        for i in range(0, ntype):
            weight_sum = np.sum(weight[i])
            if weight_sum > 0.0:
                weight[i] /= weight_sum
                tweight = np.zeros([framenum], dtype=float)
                np.copyto(tweight, weight[i])

                # calculate weighted mean and std
                ntail = len(scd[i])
                for j in range(0, ntail):
                    nc = len(scd[i][j])
                    value = np.zeros([nc, framenum], dtype=float)
                    np.copyto(value, scd[i][j])

                    for k in range(0, nc):
                        values = value[k, :]
                        ave = np.average(values, weights=tweight)
                        var = np.average((values-ave)**2, weights=tweight)
                        ascd[i][j][k], astd[i][j][k] = ave, math.sqrt(var)
                        # ascd[i][j][k],astd[i][j][k] = weighted_avg_and_std(value[k,:],tweight)
    elif otype == 'outl':
        weight = t.cast(List2D[NDFloat64], weight)
        # normalize weight
        nside = len(weight)
        for i in range(0, nside):
            ntype = len(weight[i])
            for j in range(0, ntype):
                weight_sum = np.sum(weight[i][j])
                if weight_sum > 0.0:
                    weight[i][j] /= np.sum(weight[i][j])
                    tweight = np.zeros([framenum], dtype=float)
                    np.copyto(tweight, weight[i][j])

                    # calculate weighted mean and std
                    ntail = len(scd[i][j])
                    for k in range(0, ntail):
                        nc = len(scd[i][j][k])
                        value = np.zeros([nc, framenum], dtype=float)
                        np.copyto(value, scd[i][j][k])

                        for m in range(0, nc):
                            values = value[m, :]
                            ave = np.average(values, weights=tweight)
                            var = np.average((values-ave)**2, weights=tweight)
                            ascd[i][j][k][m], astd[i][j][k][m] = ave, math.sqrt(
                                var)

    return ascd, astd
