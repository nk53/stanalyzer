#!/usr/bin/python
import sys
import typing as t

import MDAnalysis.analysis.leaflet as mdaleaflet
import numpy as np
from scipy.ndimage import gaussian_filter

if t.TYPE_CHECKING:
    import numpy.typing as npt
    from MDAnalysis import AtomGroup, Universe


NDFloat64: t.TypeAlias = "npt.NDFloat64"
nside = 2  # number of leaflets in a bilayer !!!


def find_min_dist(pos: np.ndarray, pos_grp: np.ndarray) -> 'np.floating':
    """
    ----------
    Find min dist between a position vector and a group of position vectors
    ----------

    input
      pos     : position vector
      pos_grp : position vector group (array of positions)

    output
      min_dist: minimum distance between pos and pos_grp
    """

    ngrp = len(pos_grp)  # number of vectors in the position group

    min_dist = np.linalg.norm(pos-pos_grp[0])  # initial dist
    if ngrp > 1:
        for i in range(1, ngrp):
            tdist = np.linalg.norm(pos-pos_grp[i, :])
            if tdist < min_dist:
                min_dist = tdist

    return min_dist


def assign_groups_to_leaflets(sel_grp: list['AtomGroup'],
                              ngroups: int) -> list['AtomGroup']:
    """
    ----------
    Put 3rd+ group members to the appropriate 1st or 2nd groups based on min. dist.
    ----------

    input
          sel_grp: atom group list array
          ngroups: number of atom grups in sel_grp = len(sel_grp)

    output
          sel_grp: updated sel_grp
    """

    # set arrays for job
    ng = np.zeros([ngroups], dtype=int)
    pos_grp = []  # position vector for sel_grp
    for i in range(0, ngroups):
        tmp_grp = sel_grp[i]
        ng[i] = len(tmp_grp)
        pos_grp.append(tmp_grp.positions)

    # strat from the 3rd group
    for i in range(2, ngroups):
        if ng[i] < int(0.1*ng[1]) and ng[i] < int(0.1*ng[0]):
            # Merge members into a closer group between the largest two
            for j in range(0, ng[i]):
                # find min. distances to group0 and group 1
                tpos = pos_grp[i][j, :]
                tdmin0 = find_min_dist(tpos, pos_grp[0])
                tdmin1 = find_min_dist(tpos, pos_grp[1])

                # now compare two min dist
                if tdmin0 < tdmin1:  # put this to B.groups(0)
                    # print(f'assignment {j} to group0')
                    sel_grp[0] += sel_grp[i].atoms[j]
                else:  # put this to B.groups(1)
                    # print(f'assignment {j} to group1')
                    sel_grp[1] += sel_grp[i].atoms[j]
            # print("# assignment done")
    return sel_grp


def assign_leaflet_zpos(u: 'Universe', atomgroup: 'AtomGroup') -> list['AtomGroup']:
    """
    ----------
    Assign leaflets based on z-positions of atoms.
    ----------
    Membrane normal is assumed to be parallel to the z-axis.
    zcent = (zmin + zmax)/2 is set to be a reference z-value for assignment.
    zmin and zmax are minimum and maximum z-positions of the selected atoms.

    input
          u        : MDA Universe
          atomgroup: atom group to be assinged

    output
          leaflets : atom groups for individual leaflets

    """

    leaflets: list['AtomGroup'] = []
    pos = atomgroup.positions  # position vector
    zmin = np.min(pos[:, 2])
    zmax = np.max(pos[:, 2])
    # set z value to separate leaflet
    zcent = 0.5 * (zmin + zmax)

    # upper leaflet; z > zcent
    tag0: 'AtomGroup' = atomgroup.intersection(
        u.select_atoms("prop z > %g" % zcent))
    leaflets.append(tag0)

    # lower leaflet; z < zcent
    tag1 = atomgroup.intersection(u.select_atoms("prop z < %g" % zcent))
    leaflets.append(tag1)

    # remainder - neither upper nor lower leaflet
    tag = atomgroup.difference(tag0+tag1)
    leaflets.append(tag)

    return leaflets


def assign_leaflet_mda(u: 'Universe', atomgroup: 'AtomGroup') -> list['AtomGroup']:
    """
    ----------
    Simple leaflet assignment using the default setting of LeafletFinder from MDA.
    ----------

    input
          u        : MDA Universe
          atomgroup: atom group to be assinged

    output
          leaflets : atom groups for individual leaflets

    """

    B = mdaleaflet.LeafletFinder(u, atomgroup)
    leaflets: list['AtomGroup'] = []
    for i in range(0, len(B.groups())):
        leaflets.append(B.groups(i))
        # print(f'simplest: len(B.groups({i}))= {len(leaflets[i])}')
        # if len(leaflets[i]) < 10: print(leaflets[i])

    return leaflets


def assign_leaflet_simple(u: 'Universe', atomgroup: 'AtomGroup') -> list['AtomGroup']:
    """
    ----------
    A hybrid stride-bisection method for leaflet assignment.
    ----------

    input
          u        : MDA Universe
          atomgroup: atom group to be assinged

    output
          leaflets : atom groups for individual leaflets
    """

    dcut = 15.0      # default from MDAnalysis
    dmin, dmax = (0.0, 2.0*dcut)

    # Start with default setting
    c_dcut = 0.5*(dmin+dmax)
    B = mdaleaflet.LeafletFinder(u, atomgroup, cutoff=c_dcut, pbc=True)
    ngroups = len(B.groups())

    # a hybrid stride-bisection method to obtain two groups
    if ngroups != 2:
        # initial direction
        if ngroups > 2:        # increase dmax by the width of [dmin,dmax]
            dmax += (dmax-dmin)
        # decrease dmax by the half-widht of [dmin,dmax]
        else:
            dmax -= 0.5*(dmax-dmin)

        # print(f'# iter 0, len(B.groups)={len(B.groups()} dcut={c_dcut}')

        imax = 100  # max iteration
        for i in range(1, imax+1):
            # save len(B.groups()) and cutoff distance from the previous iteration
            p_ngroups = ngroups
            p_dcut = c_dcut

            # set the current cutoff distance
            c_dcut = 0.5*(dmin+dmax)
            # assign leaflet with c_dcut
            B = mdaleaflet.LeafletFinder(u, atomgroup, cutoff=c_dcut, pbc=True)
            ngroups = len(B.groups())

            # check the current len(B.groups())
            if ngroups == 2:
                # print(f'# simple iter {i}, len(B.groups)={ngroups} after LeafletFinder, '
                #       f'dcut= {c_dcut}')
                # print("# B.groups", B.groups())
                break
            else:
                if ngroups < 2:  # decrease the cutoff distance for the next iteration
                    if p_ngroups > 2:  # p_dcut is too small
                        if p_dcut > c_dcut:
                            print(
                                f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & '
                                f'c_dcut {c_dcut} for cgr= {ngroups}')
                            sys.exit(1)
                        # update dmin,dmax
                        dmin, dmax = (p_dcut, c_dcut)
                    else:  # move along the same direction, decrease dmax
                        # decrease dmax by the half-width of [dmin,dmax]
                        dmax -= 0.5*(dmax-dmin)
                elif ngroups > 2:  # increase the cutoff distance for the next iteration
                    if p_ngroups < 2:  # p_dcut is too large
                        if p_dcut < c_dcut:
                            print(
                                f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & '
                                f'c_dcut {c_dcut} for cgr= {ngroups}')
                            sys.exit(1)
                        # update dmin,dmax
                        dmin, dmax = (c_dcut, p_dcut)
                    else:  # move along the same direction, increase dmax
                        # increase dmax by the width of [dmin,dmax]
                        dmax += (dmax-dmin)

                # print(f'# simple: iter {i}, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

        if ngroups != nside:
            print('# Job failed! Exit script !!')
            for j in range(0, len(B.groups())):
                if len(B.groups(j)) < 10:
                    print(B.groups(j))
            sys.exit(1)

    # return assigned atom groups
    leaflets = []
    for i in range(0, nside):
        leaflets.append(B.groups(i))
        # print(f'# simple: Lflt{i}: len(B.groups({i})) = {len(leaflets[i])}')

    return leaflets


def assign_leaflet(u: 'Universe', atomgroup: 'AtomGroup') -> list['AtomGroup']:
    """
    ----------
    Assign leaflet using a hybdird stride-bisection method when necessary.
    ----------
    NOTE: In typical situations that chosen atoms are appropriate for lipid tails,
          leaflets should have comparable numbers of atoms.

    input
          u        : MDA Universe
          atomgroup: atom group to be assinged

    output
          leaflets : atom groups for individual leaflets
    """

    dcut = 15.0      # default from MDAnalysis
    dmin, dmax = (0.0, 2.0*dcut)

    # Start with default setting
    c_dcut = 0.5*(dmin+dmax)
    B = mdaleaflet.LeafletFinder(u, atomgroup, cutoff=c_dcut, pbc=True)
    ngroups = len(B.groups())
    #
    sel_grp: list['AtomGroup'] = []
    ng_sel: list[int] = []
    for ig in range(0, ngroups):
        sel_grp.append(B.groups(ig))
        ng_sel.append(len(sel_grp[ig]))

    # check if selected groups for the leaflets are reasonable
    if ngroups == 2:
        if ng_sel[1] > int(0.5*ng_sel[0]):
            qdone = True
        else:
            qdone = False  # suspicious case; more groups may be assigned with smaller dcut
    elif ngroups > 2:
        # special case: the 3rd group has significantly smaller members than the 2nd group
        if ng_sel[2] < int(0.1*ng_sel[1]) and ng_sel[2] < int(0.1*ng_sel[0]):
            sel_grp = assign_groups_to_leaflets(sel_grp, ngroups)
            qdone = True
        else:
            qdone = False
    else:
        qdone = False

    # a hybrid stride-bisection method to obtain two groups
    # if ngroups != 2 and not qdone:
    if not qdone:
        # initial direction
        # suspicious case that needs decreased cutoff dist.
        if ngroups == 2:
            dmax -= 0.5*(dmax-dmin)
        elif ngroups > 2:      # increase dmax by the width of [dmin,dmax]
            dmax += (dmax-dmin)
        # decrease dmax by the half-widht of [dmin,dmax]
        else:
            dmax -= 0.5*(dmax-dmin)

        # print(f'# iter 0, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

        imax = 100  # max iteration
        for i in range(1, imax+1):
            # save len(B.groups()) and cutoff distance from the previous iteration
            p_ngroups = ngroups
            p_dcut = c_dcut

            # set the current cutoff distance
            c_dcut = 0.5*(dmin+dmax)
            # assign leaflet with c_dcut
            B = mdaleaflet.LeafletFinder(u, atomgroup, cutoff=c_dcut, pbc=True)
            ngroups = len(B.groups())
            sel_grp, ng_sel = [], []
            for ig in range(0, ngroups):
                sel_grp.append(B.groups(ig))
                ng_sel.append(len(sel_grp[ig]))

            # check the current iteration
            if ngroups == 2:
                if len(sel_grp[1]) > int(0.5*len(sel_grp[0])):
                    qdone = True
                else:
                    qdone = False

                if not qdone:
                    # suspicious case, need to decrease the cutoff distance
                    # i.e., other lipids vs. chol in the middle of flip-floping
                    # decrease dmax by the half-width of [dmin,dmax]
                    dmax -= 0.5*(dmax-dmin)
                else:
                    # print(f'# iter {i}, len(B.groups) = {ngroups} after LeafletFinder, '
                    #       f'dcut= {c_dcut}')
                    # print(f'# B.groups', B.groups())
                    break
            else:
                if ngroups > 2:  # increase the cutoff distance for the next iteration
                    # special case:
                    # the 3rd group has significantly smaller members than the 2nd group
                    if ng_sel[2] < int(0.1*ng_sel[1]) and ng_sel[2] < int(0.1*ng_sel[0]):
                        sel_grp = assign_groups_to_leaflets(sel_grp, ngroups)
                        qdone = True
                    else:
                        qdone = False
                    if qdone:
                        break

                    if p_ngroups < 2:  # p_dcut is too large
                        if p_dcut < c_dcut:
                            print(
                                f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & '
                                f'c_dcut {c_dcut} for cgr= {ngroups}')
                            sys.exit(1)
                        # update dmin,dmax
                        dmin, dmax = (c_dcut, p_dcut)
                    else:  # move along the same direction, increase dmax
                        # increase dmax by the width of [dmin,dmax]
                        dmax += (dmax-dmin)

                elif ngroups < 2:  # decrease the cutoff distance for the next iteration
                    if p_ngroups > 2:  # p_dcut is too small
                        if p_dcut > c_dcut:
                            print(
                                f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & '
                                f'c_dcut {c_dcut} for cgr= {ngroups}')
                            sys.exit(1)
                        # update dmin,dmax
                        dmin, dmax = (p_dcut, c_dcut)
                    else:  # move along the same direction, decrease dmax
                        # decrease dmax by the half-width of [dmin,dmax]
                        dmax -= 0.5*(dmax-dmin)

                # print(f'# iter {i}, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

        # if ngroups != nside:
        if qdone:
            pass
        else:
            print('# Job failed! Exit script !!')
            print(len(B.groups()))
            for j in range(0, len(B.groups())):
                if len(B.groups(j)) < 10:
                    print(B.groups(j))
                # print(B.groups(j))
            sys.exit(1)

    # print(f'# assign_leaflet: len(B.groups()) = {ngroups}')

    # return assigned atom groups
    leaflets = []
    for i in range(0, nside):
        leaflets.append(sel_grp[i])
        # print(f'# assgin_leaflet: leaflet {i}: nlipids = {len(leaflets[i])}')

    return leaflets


# ###
# #
# #  Old version: 1. atom group for atoms => 2. atom group for molecule =>
# #               3. density profile of leaflets (from molecules) =>
# #               4. find crossing point between upper & lower leaflet profiles
# #		5. return leaflets (from reprentative atoms)
# #                         z-position of the crossing point
# #
# ###
# def calc_bilayer_midplane_old(u,ag_sys,nside,method):
#   """
#   ----------
#   Calculate bilayer center based on leaflet density profiles (e.g. termnial methyl)
#   ----------
#
#   input
#         u       : MDA Universe
#         ag_sys  : atom group for membrane. Recommended to select all lipid types.
#         nside   : number of leaflets (=2)
#   output
#         leaflets: leaflet atom groups from ag_sys
#         zcent   : bilayer center along z
#
#   NOTE: currently, only PSF based analysis is supported.
#   """
#
#   # Generate leaflet atom groups
#   if method == 'zpos':
#     leaflets = assign_leaflet_zpos(u,ag_sys)
#   elif method == 'mda':
#     leaflets = assign_leaflet(u,ag_sys)
#
#   # Process to get heavy atoms of membrane
#   ag_memb=[]
#   for i in range(0,nside):
#     ag_memb.append(u.atoms[[]]) # empty atom group
#     natom = len(leaflets[i])    # number of atoms
#     seg = leaflets[i].segids
#     rid = leaflets[i].resids
#     rnm = leaflets[i].resnames
#     # new selections & atom groups
#     for j in range(0,natom):
#       tag = u.select_atoms(f'segid {seg[j]} and resname {rnm[j]} and resid {rid[j]} '
#                            'and not name H*')
#       ag_memb[i] += tag
#
#   # Calculate leaflet density profiles
#   pos=[]
#   for i in range(0,nside): pos.append(ag_memb[i].positions)
#
#   # get histogram of heavy atom number z-profile
#   zmin,zmax,dz=(-20.0,20.0,0.4) # dz=0.1 is too fine...
#   nbin=int(2.0*zmax/dz)
#   zbins=[]
#   for i in range(0,nbin+1): zbins.append(zmin+i*dz)
#
#   hist = []
#   for i in range(0,nside):
#     hist.append([])
#     hist[i],bin_edges = np.histogram(pos[i][:,2],bins=zbins)
#
#   # Smoothing histogram
#   sig = int(1.0/dz) # Set sig = 1.0 A
#   if sig < 1.0: sig = 1.0
#   for i in range(0,nside):
#     hist[i]=gaussian_filter(hist[i],sigma=sig)
#
#   # Find index where Xing occurs
#   idx = np.argwhere(np.diff(np.sign(hist[0]-hist[1]))).flatten()
#
#   # Estimate zcenter from the intersection of two line segments
#   i0,i1 = idx[0], int(idx[0]+1)
#   fx0,fx1 = (hist[0][i0],hist[0][i1])
#   gx0,gx1 = (hist[1][i0],hist[1][i1])
#   zcent = zmin + (i0+0.5)*dz - (fx0-gx0)*dz/((fx1-fx0)-(gx1-gx0))
#
#   return leaflets,zcent


def calc_bilayer_midplane(ag_memb, nside):
    """
    ----------
    Calculate bilayer center based on leaflet density profiles (e.g. all heavy atoms)
    ----------

    input
          ag_memb : atom groups for membrane. Assummed to be lipid tails or molecules
          nside   : number of leaflets (=2)
    output
          zcent   : bilayer center along z
    """

    # Calculate leaflet density profiles
    pos = []
    for i in range(0, nside):
        pos.append(ag_memb[i].positions)

    # get histogram of heavy atom number z-profile
    zmin, zmax, dz = (-20.0, 20.0, 0.4)  # dz=0.1 is too fine...
    nbin = int(2.0*zmax/dz)
    zbins = []
    for i in range(0, nbin+1):
        zbins.append(zmin+i*dz)

    hist: list[NDFloat64] = []
    for i in range(0, nside):
        hist.append([])
        hist[i], bin_edges = np.histogram(pos[i][:, 2], bins=zbins)

    # Smoothing histogram
    sig = float(int(1.0/dz))  # Set sig = 1.0 A
    if sig < 1.0:
        sig = 1.0
    for i in range(0, nside):
        hist[i] = gaussian_filter(hist[i], sigma=sig)

    # Find index where Xing occurs
    idx = np.argwhere(np.diff(np.sign(hist[0]-hist[1]))).flatten()

    # Estimate zcenter from the intersection of two line segments
    i0, i1 = idx[0], int(idx[0]+1)
    fx0, fx1 = (hist[0][i0], hist[0][i1])
    gx0, gx1 = (hist[1][i0], hist[1][i1])
    zcent = zmin + (i0+0.5)*dz - (fx0-gx0)*dz/((fx1-fx0)-(gx1-gx0))

    return zcent
