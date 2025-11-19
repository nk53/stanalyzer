#!/usr/bin/python
import sys
import typing as t

import numpy as np
from . import leaflet_util as myleaflet

if t.TYPE_CHECKING:
    from MDAnalysis import AtomGroup, Universe


class MolTypeAGs(t.NamedTuple):
    sel_name: str
    sel_ag: list['AtomGroup']
    sel_nmol: int


class MolGroups(t.NamedTuple):
    name_type: list[str]
    nmol_type: list[int]
    nmol: int
    id_type: list[int]
    ag: list['AtomGroup']


class MolGroupsMemb(t.NamedTuple):
    name_type: list[str]
    nmol_type: list[list[int]]
    nmol: list[int]
    id_type: list[list[int]]
    ag: list['AtomGroup']


class MolGroupsFull(t.NamedTuple):
    nmol_type: list[int]
    nmol: int
    id_type: list[int]
    ag_full: list['AtomGroup']


def generate_mol_type_ags(ag: 'AtomGroup', selection: str, qsplit: bool) -> MolTypeAGs:
    """
    ----------
    Generate atom groups for the given selection of a molecule type
    ----------
    NOTE: Get atom group & split into molecules based on selection & qsplit.

    input
          ag       : atomgroup
          selection: atom selection for a given molecule type
          qsplit   : True (split to molecule level)/ False (treat selection as a molecule)

    output
          nmol     : number of molecules for the selection
          sel_ag   : atom groups for individual molecules
    """

    tmps = selection.split()
    # selection starts with either segid, resname, or moleculetype
    # moleculetype is a specific feature for GROMACS
    # u.atoms.select("moleculetype 'moleculename'")
    #  -> molnums => molecule id, mol type => molecule names

    sel_type = tmps[0]  # segid/resname/moleculetype
    sel_name = tmps[1]  # PRO*/DSPC/DOPC/ ...

    if sel_type == "segid":
        if qsplit:
            sel_ag = ag.split('segment')
        else:  # if (qsplit is False): # no split
            sel_ag = []
            sel_ag.append(ag)
        # for i in range(0,len(sel_ag)):
        #   print(sel_ag[i][0].segid,sel_ag[i][0].resid,sel_ag[i][0].resname)
    elif sel_type == "resname":
        if qsplit:
            sel_ag = ag.split('residue')
        else:
            sel_ag = []
            sel_ag.append(ag)
    # NEED to test!!!
    elif sel_type == "moltype":  # GROMACS specific
        if qsplit:
            sel_ag = ag.split('molecule')
        else:
            sel_ag = []
            sel_ag.append(ag.molecules)
    else:  # neither segid nor resname
        print('# unique molecule type is not defined as moleculetype/segid/resname')
        sys.exit(1)

    sel_nmol = len(sel_ag)

    return MolTypeAGs(sel_name, sel_ag, sel_nmol)


def generate_mol_groups(u: 'Universe', ntype: int, selection: list[str],
                        qsplit: list[bool]) -> MolGroups:
    """
    ----------
    Generate atom groups for individual molecules.
    ----------

    input
          u         : MDA Universe
          ntype     : number of molecule types
          selection : selections for individual molecule types
          qsplit    : option for splitting atom groups into molecule level

    output
          name_types: names of individual molecule types
          nmol_type : numbers of molecules for individual molecule types
          nmol      : total number of molecules
          id_type   : molecule type indices for individual molecules
          ag        : atom groups for individual molecules
    """

    name_type: list[str] = []  # names of unique molecule types
    nmol_type: list[int] = []  # numbers of individual molecule types
    id_type: list[int] = []  # molecule type indices of individual molecules

    # generate atom groups for individual molecules
    imol = 0
    ag: list['AtomGroup'] = []
    for i in range(0, ntype):
        ag_tmp = u.select_atoms(selection[i])

        nt, sel_ag, sel_nmol = generate_mol_type_ags(
            ag_tmp, selection[i], qsplit[i])
        name_type.append(nt)

        nmol_type.append(sel_nmol)  # number of molecules of type, i
        for j in range(0, sel_nmol):
            # append atom group for individual molecule
            ag.append(u.atoms[[]])
            ag[imol] += sel_ag[j]     # update molecular atom group
            id_type.append(i)         # update molecule type index
            imol += 1                 # update molecule index counter

    nmol: int = np.sum(nmol_type)      # total number of molecules

    # print system info
    sout = '# SYS. INFO\n'
    sout += '# number of molecule types = {ntype}\n'
    sout += '# type  :'
    for i in range(0, ntype):
        sout += f' {name_type[i]:10s}'
    sout += '\n'
    sout += '# index :'
    for i in range(0, ntype):
        sout += f' {i:10d}'
    sout += '\n'
    sout += '# number:'
    for i in range(0, ntype):
        sout += f' {nmol_type[i]:10d}'
    sout += f' # nmol = {nmol:5d} ; # of ags ={len(ag)}\n'
    print(sout)

    return MolGroups(name_type, nmol_type, nmol, id_type, ag)


def generate_mol_groups_memb(u: 'Universe', nside: int, ntype: int,
                             selection: list[str], qsplit: list[bool], sside: list[str],
                             qz: bool) -> MolGroupsMemb:
    """
    ----------
    Generate atom groups for individual molecules in individual leaflets.
    ----------
    input
          u         : MDA Universe
          nside     : number of leaflets (=2)
          ntype     : number of molecule types
          selection : selections for individual molecule types
          qsplit    : split option for individual selections
          sside     : leaflet names (["up","dn"]).
          qz        : option for leaflet assignemnt.
                      If True, z-position based assignment, otherwise MDA leaflet.

    output
          name_types: names of individual molecule types
          nmol_type : numbers of molecules for individual molecule types
          nmol      : total number of molecules in individual leaflets
          id_type   : molecule type indices for individual molecules in indv. leaflets
          ag        : atom groups for individual molecules in individual leaflets
    """

    name_type: list[str] = []  # names of unique molecule types
    nmol_type: list[list[int]] = []  # numbers of individual molecule types
    nmol: list[int] = []  # number of molecules
    id_type: list[list[int]] = []  # molecule type indices of individual molecules

    # expand nmol_type and id_type for leaflets
    for i in range(0, nside):
        nmol_type.append([])
        nmol.append(-1)
        id_type.append([])

    # generate ordered atom groups for individual lipid types
    # an empty atom group to which ag_type[] will be added
    ag_tmp = u.atoms[[]]
    for i in range(ntype):
        ag_tmp += u.select_atoms(selection[i])

    # Assign leaflet at the beginning & use this for MSD
    # Do not consider flip-flop! - It cannot be handled properly
    if qz:
        print('# leaflet assignment based on z-position')
        leaflet = myleaflet.assign_leaflet_zpos(u, ag_tmp)
    else:
        print('# leaflet assignemnt based on LeafletFinder with hydrid cutoff search')
        leaflet = myleaflet.assign_leaflet(u, ag_tmp)
    # print(leaflet)
    # sys.exit(0)

    # now setup leaflet ag
    ag: list['AtomGroup'] = []
    for i in range(0, nside):
        ag.append([])
        imol = 0

        # loop over ntypes
        for j in range(0, ntype):
            # atom group for selection in leaflet[i]
            ag_tmp = u.select_atoms(selection[j]).intersection(leaflet[i])

            # split atom group into segid/resname/moleculetype (i.e., molecule level)
            tname_type, sel_ag, sel_nmol = generate_mol_type_ags(
                ag_tmp, selection[j], qsplit[j])

            if i == 0:
                name_type.append(tname_type)

            nmol_type[i].append(sel_nmol)  # number of molecules of type, i
            for k in range(0, sel_nmol):
                # append atomgroup for individual molecule
                ag[i].append(u.atoms[[]])
                ag[i][imol] += sel_ag[k]     # update molecular atom group
                id_type[i].append(j)         # update molecule type index
                imol += 1                    # update molecule index counter

        # total number of molecules in the leaflet, i
        nmol[i] = np.sum(nmol_type[i])

    # print system info
    sout = '# SYS. INFO\n'
    sout += f'# number of molecule types = {ntype}\n'
    sout += '# type  :'
    for i in range(0, ntype):
        sout += f' {name_type[i]:10s}'
    sout += '\n'
    sout += '# index :'
    for i in range(0, ntype):
        sout += f' {i:10d}'
    sout += '\n'
    sout += '# number\n'

    for i in range(0, nside):
        sout += f'# {sside[i]:7s}'
        for j in range(0, ntype):
            sout += f' {nmol_type[i][j]:10d}'
        sout += f' # nmol = {nmol[i]:5d} ; # of ags = {len(ag[i]):5d}\n'
    print(sout)

    return MolGroupsMemb(name_type, nmol_type, nmol, id_type, ag)


def generate_full_mol_groups(u: 'Universe', ntype: int, sel_type: list[str],
                             name_type: list[str], qsplit: list[bool]) -> MolGroupsFull:
    """
    ----------
    Generate reference atomgroups for individual molecules with full atoms.
    ----------
    NOTE: It can be used to molecule identification by intersectioning other atom groups.

    input
          u        : MDA Universe
          ntype    : numbers of molecule types
          sel_type : selection type; segid/resname/moltype
          name_type: names of individual molecule types
          qsplit   : options for splitting selections into molecule level

    output
          nmol_type: numbers of molecules for individual molecule types
          nmol     : total number of molecules
          id_type  : molecule type indices for individaul molecules
          ag_full  : atom groups of individual molecules with full atoms
    """

    nmol_type = []  # numbers of individual molecule types
    id_type = []  # molecule type indices of individual molecules

    # generate atom groups for individual molecules
    imol, ag_full = 0, []
    for i in range(0, ntype):
        selection = f'{sel_type[i]} {name_type[i]}'
        ag_tmp = u.select_atoms(selection)  # single atom selection

        tname_type, sel_ag, sel_nmol = generate_mol_type_ags(
            ag_tmp, selection, qsplit[i])

        nmol_type.append(sel_nmol)  # number of molecules of type, i
        for j in range(0, sel_nmol):
            # append atom group for individual molecule
            ag_full.append(u.atoms[[]])
            ag_full[imol] += sel_ag[j]     # update molecular atom group
            id_type.append(i)         # update molecule type index
            imol += 1                 # update molecule index counter

    nmol = np.sum(nmol_type)      # total number of molecules

    # print system info
    sout = '# SYS. INFO\n'
    sout += f'# number of molecule types = {ntype}\n'
    sout += '# type  :'
    for i in range(0, ntype):
        sout += f' {name_type[i]:10s}'
    sout += '\n'
    sout += '# index :'
    for i in range(0, ntype):
        sout += f' {i:10d}'
    sout += '\n'
    sout += '# number:'
    for i in range(0, ntype):
        sout += f' {nmol_type[i]:10d}'
    sout += f' # nmol = {nmol:5d} ; # of ags = {len(ag_full)}\n'
    print(sout)

    return MolGroupsFull(nmol_type, nmol, id_type, ag_full)


def generate_sys_groups(u: 'Universe', sel_sys: str, qz: bool) -> list['AtomGroup']:
    """
    ----------
    Generate system atom groups for individual leaflets
    ----------
    NOTE: Used for the leaflet COM drift correction in MSD analysis of lipids.

    input
          u      : MDA Universe
          sel_sys: selection for system
          qz     : option for leaflet assignemnt.
                   If True, z-position based assignment, otherwise use MDA leaflet.

    output
          ag_sys : atomgroups of system for individual leaflets
    """

    ag_tmp = u.select_atoms(sel_sys)
    # Leaflet assignment
    if qz:
        ag_sys = myleaflet.assign_leaflet_zpos(u, ag_tmp)
    else:
        ag_sys = myleaflet.assign_leaflet(u, ag_tmp)

    return ag_sys
