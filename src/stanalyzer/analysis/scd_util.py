#!/usr/bin/python
import sys
import math
import re
import numpy as np

# nside=2     # up/dn


def find_carbons_hydrogens(u, atomgroup, name_t, c_st, ntail):
    """
    ----------
    Find carbonds and associated bonded hydrogens in tails of a given molecule
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
    carbons = []
    hydrogens = []
    ncarbons = []

    # expand arrays along chains
    for i in range(0, ntail):
        carbons.append([])
        hydrogens.append([])
        ncarbons.append(0)

    # Search for carbons
    for i in range(0, ntail):
        c_curr = c_st[i]  # starting carbon name in the chain

        # start from ic_name of each chain until there exist no more connected carbon atoms to the c_curr
        tag_c = atomgroup.intersection(u.select_atoms(f'name {c_curr}'))
        tag_h = atomgroup.intersection(
            u.select_atoms(f'name H* and bonded name {c_curr}'))

        if len(tag_c) == 0:
            print(f'# Starging Carbon name {c_curr} does not exist!')
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
                print(f'# {name_t} chain {i}: all carbon atoms assigned')
                break

        # Print out assigned carbon/bonded hydrogens
        for j in range(0, ncarbons[i]):
            print(f'# {name_t}: chain {i} ica = {j}',
                  carbons[i][j], hydrogens[i][j])

    return carbons, hydrogens, ncarbons


def generate_carbons_hydrogens_ncarbons_type(u, ag_full, ntype, name_type,
                                             nmol_type, ic_name, nchain):
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

    carbons = []    # [ntype][nchain][0-ncarbons]
    hydrogens = []  # [ntype][nchain][0-ncarbons][...]
    ncarbons = []

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

    return carbons, hydrogens, ncarbons


def generate_ag_carbon_hydrogen(u, atomgroups, name_type, id_type, nchain,
                                ncarbons, carbons, hydrogens):
    """
    ----------
    Generate atom groups of carbon/associated bonded hydrogens for individual lipids
    ----------

    input
          u         : MDA Universe
          atomgroups: atomgroups of individual molecules
          nchain    : number of chains in individual molecule types
          carbons   : carbon names in ind. chain in unique molecule types
          hydrogens : hydrogen names bonded to ind. carbons in ind. chains
          carbons   : number of carbons in ind. chains in unique mol. types

    output
          ag_c      : atom groups of carbons in ind. chain in ind. lipids
          ag_h      : atom groups of Hs bonded to ind. Cs in ind. tail in ind. lipids
    """

    nmol = len(atomgroups)

    # atom group list arrays
    ag_c, ag_h = [], []

    # Loop over lidids
    for i in range(0, nmol):
        ag_c.append([])
        ag_h.append([])

        itype = id_type[i]        # lipid type index
        name_t = name_type[itype]  # molecule type name

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


def generate_ref_atomgroups(u, atomgroups, ntype, nmol_type, ref_name):
    """
    ----------
    Generate (ordered) reference atom groups for leaflet assignment
    ----------

    input
          u         : MDA Universe
          atomgroups: atom groups of all molecules with all atoms
          ntype     : number of molecule types
          nmol_type : numbers of molecules of individual molecule types
          ref_name  : reference atom name for individual molecule types

    output
          ag_ref    : reference atom groups
    """

    ag_ref = u.atoms[[]]
    imol = 0
    for i in range(0, ntype):
        for j in range(0, nmol_type[i]):
            sel_ref = f'name {ref_name[i]}'
            tag = atomgroups[imol].intersection(u.select_atoms(sel_ref))
            # ag_ref.append(tag)
            ag_ref += tag
            imol += 1
    # print(ag_ref)
    return ag_ref


def setup_raw_scd(nmol, id_type, nchain, ncarbons):
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

    raw_scd = []
    for i in range(0, nmol):
        # lipid type and associated number of tails
        itype = id_type[i]
        ntail = nchain[itype]

        raw_scd.append([])
        for j in range(0, ntail):
            nc = ncarbons[itype][j]

            # now make numpy array
            raw_scd[i].append(np.zeros([nc], dtype=float))
    return raw_scd


def calculate_raw_scd(raw_scd, nmol, ag_c, ag_h, memb_norm):
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


def setup_scd_weight(nside, ntype, nchain, ncarbons, framenum, qb):
    """
    ----------
    Set up SCD time series arrays and associated un-normalized weights
    ----------

    input
          nside   : number of leaflet
          ntype   : number of unique lipid types
          nchain  : number of chains in unique lipid types
          ncarbons: number of carbons in ind. chain in unique lipid types
          qb      : False (default: leaflet analysis)/True (bilayer analysis)

    output
          scd     : SCD of ind. carbon in ind. chain in unique lipid types (bilayer/leaflet)
          weight  : number of counts for unique lipid types (bilayer/leaflet)
    """

    scd = []
    weight = []  # number of counts
    if qb:
        for i in range(0, ntype):
            scd.append([])
            weight.append(np.zeros([framenum], dtype=float))

            ntail = nchain[i]
            for j in range(0, ntail):
                nc = ncarbons[i][j]
                scd[i].append(np.zeros([nc, framenum], dtype=float))
    else:  # if (qb is False):

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

    print(np.shape(scd))
    return scd, weight


def setup_ascd_astd(nside, ntype, nchain, ncarbons, qb):
    """
    ----------
    Setup average and std of SCDs
    ----------

    input
          nside   : number of leaflets
          ntype   : number of unique lipid types
          nchain  : number of chains in unique lipid types
          ncarbons: number of carbons in ind. chains in ind. lipid types
          qb      : True (bilayer analysis)/ False (default: leaflet analysis)

    output
          ascd    : average SCD
          astd    : std of SCD
    """

    ascd = []  # weighted mean
    astd = []  # weighted std
    if qb:
        for i in range(0, ntype):
            ascd.append([])
            astd.append([])

            ntail = nchain[i]
            for j in range(0, ntail):
                nc = ncarbons[i][j]
                ascd[i].append(np.zeros([nc], dtype=float))
                astd[i].append(np.zeros([nc], dtype=float))
    else:  # if (qb is False):
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

    print(np.shape(ascd))
    return ascd, astd


def assign_leaflet_index(ag_ref, leaflets):
    """
    ----------
    Assign leaflet index to individual lipids (for a given frame)
    ----------

    input
          ag_ref  : reference atomgroup for leaflet assignment
          leaflets: list array of individual leaflets from ag_ref

    output
          id_side : list array of leaflet index of individual lipids
    """

    nmol = len(ag_ref)

    tid_side = []
    [tid_side.append(-1) for i in range(0, nmol)]
    flag = []
    [flag.append(0) for i in range(0, nmol)]  # before assignment
    nside = len(leaflets)  # number of leaflets

    # Loop over molecule
    for i in range(0, nmol):
        for j in range(0, nside):
            if flag[i] == 0:
                tag = leaflets[j].intersection(ag_ref[i])
                if len(tag) > 0:
                    # print(tag)
                    tid_side[i] = j
                    flag[i] = 1
                    break

    # check for unassigned lipid
    for i in range(0, nmol):
        if flag[i] == 0:
            print(f'# unassigned lipid, {i}')
            sys.exit(1)

    return tid_side


def calculate_scd(iframe, ntype, raw_scd, scd, weight, id_side, id_type, nchain, ncarbons, qb):
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
          qb      : True (bilayer)/ False (leaflet)

    output
          scd,weight
    """

    nmol = len(raw_scd)
    # ntype: global parameter
    # nsdie: global parameter

    if qb:
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

    else:  # if (qb is False):
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


def calculate_ave_and_std(ascd, astd, scd, weight, framenum, qb):
    """
    ----------
    Calculate average and std of SCDs for individaul lipid types
    ----------
    input
          scd   : time series of SCD for ind. carbons in ind. chains in unique lipid types
          weight: un-normalized wight for each time frame
          qb    : True (bilayer)/ False (Leaflet)

    input/output
          ascd  : weighted ave SCD for ind. carbons in ind. chains in unique lipid types
          astd  : weighted std of SCD
    """

    # nside/ntype : global parameters
    if qb:
        # normalize weight
        ntype = len(weight)
        for i in range(0, ntype):
            weight[i] /= np.sum(weight[i])
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
    else:  # if (qb is False):
        # normalize weight
        nside = len(weight)
        for i in range(0, nside):
            ntype = len(weight[i])
            for j in range(0, ntype):
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
