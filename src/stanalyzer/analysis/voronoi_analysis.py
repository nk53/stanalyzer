#!/usr/bin/python
import sys
import math
import typing as t
from collections.abc import Sequence

import numpy as np

if t.TYPE_CHECKING:
    import numpy.typing as npt
    from scipy.spatial import Voronoi
    from MDAnalysis import AtomGroup, Universe


T = t.TypeVar('T')
List2D = list[list[T]]
List3D = list[list[list[T]]]
List4D = list[list[list[list[T]]]]
List5D = list[list[list[list[list[T]]]]]
NDFloat64: t.TypeAlias = 'npt.NDArray[np.float64]'
NDInt64: t.TypeAlias = 'npt.NDArray[np.int64]'

OutputFileType: t.TypeAlias = t.Literal['outl', 'outb']

# --- The following are hard set for voronoi tessellation in 2D
# dim=2              # 2 dimensions
# nimage=int(3**dim) # number of total primary+images for 2D


def connect_mol_and_vc(u: 'Universe', nmol: int, nmol_type: list[int],
                       ag: list['AtomGroup']) -> tuple[List2D[int], NDInt64]:
    """
    ----------
    Generate VC list for molecules &
    assign molecule index to individual VCs
    ----------

    input
          u        : MDA Universe
          nmol     : number of molecules
          nmol_type: numbers of molecules of individual molecule types
          ag       : atom groups of individual molecules

    output
          mol_ivc  : list of VCs for individual molecules
          id_mol_vc: molecule index for individual VCs
    """

    # total number of VCs in the atom groups, ag
    npt = 0
    for i in range(0, nmol):
        npt += len(ag[i])
    print(f'# total number of VCs = {npt}')

    id_mol_vc = np.zeros([npt], dtype=int)
    mol_ivc: List2D[int] = []

    iatom = 0
    for i in range(0, nmol):
        mol_ivc.append([])
        natom = len(ag[i])
        for j in range(0, natom):
            mol_ivc[i].append(iatom)
            id_mol_vc[iatom] = i
            iatom += 1

    return mol_ivc, id_mol_vc


def image_translation(box: NDFloat64) -> NDFloat64:
    """
    ----------
    Translation operators for images
    ----------

    input
          box : XY-box dimension vector

    output
          tran: array of translation vector for individual images
    """

    # box dimensions
    xtla, xtlb = box[0], box[1]  # noqa: F841

    # translation sequence (primary, ..., last image)
    # Hard set operator for image translation in 2D
    operator = np.array([[0, 0],  [1, 0],  [1, 1],
                         [0, 1],  [-1, 1], [-1, 0],
                         [-1, -1], [0, -1], [1, -1]], dtype=float)
    noperator = len(operator)  # ! should be equal to nimage

    # translation vectors
    tran = np.zeros([noperator, 2], dtype=float)
    for i in range(0, noperator):
        np.copyto(tran[i, :], operator[i, :]*box)

    return tran


def prepare_voronoi_centers(ag: list['AtomGroup'], box: NDFloat64, dim: int,
                            nimage: int) -> tuple[int, int, NDFloat64]:
    """
    ----------
    Prepare voronoi centers (primary + 8 images)
    ----------

    input
          ag    : atom groups for molecules in a chosen leaflet
          box   : xy-box dim vector
          dim   : 2
          nimage:

    output
          nprim: number of primary Voronoi centers
          nvorn: number of total Voronoi centers
          vcent: voronoi centers
    """

    # set number of VCs for primary simulation box points
    nprim = 0
    nmol = len(ag)  # atom groups for individual molecules
    for i in range(0, nmol):
        nprim += len(ag[i])
    nvorn = int(nimage*nprim)    # total VCs in 2D (primary + images)

    sout = f'# DIMENSION= {dim} NIMAGE= {nimage}\n'
    sout += '# VORONOI CENTER INFO.\n'
    sout += f'       nprim= {nprim} nvorn= {nvorn}\n'
    print(sout)

    # initialize VC array
    vcent = np.zeros([nimage*nprim, dim], dtype=float)
    tvcen = np.zeros([nprim, dim], dtype=float)

    iatom = 0
    for i in range(0, nmol):
        natom = len(ag[i])  # number of atoms
        ist, ied = (iatom, iatom+natom)  # start/end indices in tvcen
        tpos = ag[i].positions
        # copy xy coordinates to tvcen
        np.copyto(tvcen[ist:ied], tpos[:, 0:dim])
        iatom += natom

    # translation vector for all images
    tran = image_translation(box)

    # loop over images
    for iimg in range(0, nimage):

        tmpvc = tvcen + tran[iimg]  # translation

        # start & stop indices in vcent array
        ist, ied = int(iimg*nprim), int((iimg+1)*nprim)

        # copy image coordinates to vcent array
        np.copyto(vcent[ist:ied], tmpvc)
        # print(f'# iimg={iimg}, xtla={box[0]} xtlb={box[1]} vcent[ist]= {vcent[ist]}')

    return nprim, nvorn, vcent


# --------------------
# Functions for CONTACT analysis
# --------------------


def generate_vc_contact_info(vorn: 'Voronoi', nprim: int, nimage: int, box: NDFloat64) \
        -> tuple[list[NDInt64], list[NDFloat64], list[NDFloat64], list[NDFloat64]]:
    """
    ----------
    Generate a list array of contacting VCs for individual VCs in the primary cell
    ----------

    input
          vorn     : Voronoi(atomgroup)
          nprim    : number of atoms in the primary image
          nimage   : number of images, set to be 9 (for 2D)
          box      : box vector

    output
          contactVC: a list of ndarray of the list of contacting VCs to individual VC
          borderVC : a list of ndarray of the border btw. individual VC and its NN VCs
          distVC   : a list of ndarray of the distances btw. individual VC and its NN VCs
          weightVC : a list of ndarray of the weight=dist/border btw. ind. VC and its NN VCs
    """

    # set contactVC for the VCs in the primary image
    _contactVC: list[list[int]] = [[] for _ in range(nprim)]
    _borderVC:  list[list[np.floating]] = [[] for _ in range(nprim)]
    _distVC:    list[list[np.floating]] = [[] for _ in range(nprim)]
    _weightVC:  list[list[np.floating]] = [[] for _ in range(nprim)]

    # search for vorn.ridge_points
    nridges = len(vorn.ridge_points)
    for i in range(0, nridges):
        inx1: int
        inx2: int
        inx1, inx2 = vorn.ridge_points[i, 0], vorn.ridge_points[i, 1]

        # check validity
        if inx1 < 0 or inx2 < 0:
            print(f'# construct N.N. VC: Error - index range: {inx1} {inx2}')
            sys.exit(1)
        # generate only for the VCs in the primary image
        if inx1 >= nprim and inx2 >= nprim:
            continue

        # inx1 < nprim or inx2 < nprim
        # calc border length btw. VCs
        # indices of the voronoi vertices for the ridge
        vinx = vorn.ridge_vertices[i]
        vert1 = vorn.vertices[vinx[0]]
        vert2 = vorn.vertices[vinx[1]]
        tborder = np.linalg.norm(vert1 - vert2)
        # calc dist btw. VCs
        tdist = np.linalg.norm(vorn.points[inx1]-vorn.points[inx2])

        # 1.  add indices to contactVC & update distVC and borderVC
        if inx1 < nprim:
            _contactVC[inx1].append(inx2)
            _borderVC[inx1].append(tborder)
            _distVC[inx1].append(tdist)
            _weightVC[inx1].append(tborder/tdist)

        if inx2 < nprim:
            _contactVC[inx2].append(inx1)
            _borderVC[inx2].append(tborder)
            _distVC[inx2].append(tdist)
            _weightVC[inx2].append(tborder/tdist)

    # convert into numpy array
    contactVC: list[NDInt64] = [np.array(tarray) for tarray in _contactVC]
    borderVC:  list[NDFloat64] = [np.array(tarray) for tarray in _borderVC]
    distVC:    list[NDFloat64] = [np.array(tarray) for tarray in _distVC]
    weightVC:  list[NDFloat64] = [np.array(tarray) for tarray in _weightVC]

    return contactVC, borderVC, distVC, weightVC


def generate_mol_contact_info(mol_ivc: list[list[int]],
                              contactVC: list[NDInt64],
                              id_mol_vc: NDInt64,
                              nprim: int) -> list[set[int]]:
    """
    ----------
    Generate a list array of molecular contact for individual molecules in the primary cell
    ----------

    input
          id_mol_vc: molecule index for individual VC
          mol_ivc  : VC list for individual molecule
          contactVC: neighbor VC list for individual VC
          nprim    : number of VCs in the primary image

    output
          contactMol: a list array of set of neighboring molecules for individual molecules
    """

    # number of molecules
    nmol = len(mol_ivc)

    # empyt list output array
    contactMol = []

    # Loop over molecules
    for i in range(0, nmol):
        # collect all contactVCs for molecule, i
        tset = set()  # empty set for the contactVCs for molecule, i
        nsite = len(mol_ivc[i])  # noqa: F841
        for elem1 in mol_ivc[i]:
            for elem2 in contactVC[elem1]:
                tid: int = id_mol_vc[elem2 % nprim]  # molecule index
                if tid != i:
                    tset.add(tid)

        # append to the output array
        contactMol.append(tset)

    return contactMol


def calculate_contact(
        nmol: int,
        ntype: int,
        nmol_type: list[int],
        mol_ivc: List2D[int],
        id_type: list[int],
        contactMol: list[set[int]]) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Calculate contacting lipid types for individual lipid types
    ----------

    input
          nmol      : number of molecules
          ntype     : number of lipid types
          nmol_type : numbers of individual lipid types
          mol_ivc   : list of VCs for individual molecules
          id_type   : molecule type indices for individual molecules
          contactMol: list of contacting molecules for individual molecules

    output
          ncontact  : numbers of contacting mol. types for ind. mol. types
          fcontact  : fractions of contacting mol. types for ind. mol. typ
    """

    # define ncontact and fcontact
    ncontact: NDFloat64 = np.zeros([ntype, ntype], dtype=float)
    fcontact: NDFloat64 = np.zeros([ntype, ntype], dtype=float)

    # loop over mols
    for i in range(0, nmol):
        inx1 = id_type[i]
        tnn_mol = contactMol[i]
        for j in tnn_mol:
            inx2 = id_type[j]
            ncontact[inx1][inx2] += 1.0

    # division by the number of components
    for i in range(0, ntype):
        norm = nmol_type[i]
        if norm > 0.0:
            ncontact[i, :] /= norm

    # calculate fcontact
    for i in range(0, ntype):
        tnorm = np.sum(ncontact[i, :])
        if tnorm > 0.0:
            fcontact[i] = ncontact[i]/tnorm

    return ncontact, fcontact

# --------------
#
# Update contact for bilayer from leaflet data
#
# --------------


def update_contact_bilayer(framenum: int, ntype: int, ncontact: NDFloat64, fcontact: NDFloat64,
                           tnmol_type: NDFloat64) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Update contacts for bilayer: weighted average of leaflet data
    ----------

    input
          framenum  : total number of frames in input trajectories
          ntype     : number of molecule types
          tnmol_type: nomrmalization factor (sum of weights over leaflets)

    input/output
          ncontact  : time series of contact numbers
          fcontact  : time series of fractional contacts
    """

    for i in range(0, framenum):
        # update ncontact
        for j in range(0, ntype):
            tnorm = tnmol_type[i, j]  # it's component should not be zero.
            ncontact[i, j, :] /= tnorm
        # update fcontact
        for j in range(0, ntype):
            tnorm = np.sum(ncontact[i, j, :])
            if tnorm > 0.0:
                fcontact[i, j, :] = ncontact[i, j, :]/tnorm

    return ncontact, fcontact


# --------------------
# Functions for Shell analysis
# --------------------

def setup_shell_comp_arrays(framenum: int, nside: int, ntype: int,
                            max_shell: int, otype: OutputFileType) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Setup raw shell-wise composition arrays
    ----------

    input
          framenum : number of frames for the analysis
          nside    : 2 (bilayer)
          ntype    : number of unique molecule types
          otype    : Output type: outl (Leaflets); outb (Bilayer)

    output
          t_numb_comp: list array of time series of shell-wise number compositions
          t_frac_comp: list array of time series of shell-wise fractional compositions
    """

    def do_small() -> tuple[NDFloat64, NDFloat64]:
        t_numb_comp = np.zeros(
            [framenum, ntype, max_shell+1, ntype], dtype=float)
        t_frac_comp = np.zeros(
            [framenum, ntype, max_shell+1, ntype], dtype=float)
        return t_numb_comp, t_frac_comp

    def do_big() -> tuple[NDFloat64, NDFloat64]:
        t_numb_comp = np.zeros(
            [framenum, nside, ntype, max_shell+1, ntype], dtype=float)
        t_frac_comp = np.zeros(
            [framenum, nside, ntype, max_shell+1, ntype], dtype=float)
        return t_numb_comp, t_frac_comp

    if otype == 'outb':
        return do_small()
    return do_big()


def calculate_shell_comp(imol: int,
                         nmol: int,
                         contactMol: list[set[int]],
                         ntype: int,
                         mol_ivc: list[list[int]],
                         id_type: list[int],
                         max_shell: int) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Calculate shell-wise number and fractional compositions around a molecule
    up to max_shell-th shell (including images)
    ----------

    input
          imol      : molecule index to generate shell arount it
          nmol      : total number of molecules
          contactMol: contacting molecule index list for individual molecule
          mol_ivc   : voronoi cell list for individual molecules
          id_type   : lipid type index for individual molecule
          max_shell : max shell index

    output
          ARRAYS of [max_shell,ntype]: NDFloat64
          m_numb_comp : numbers of lipid types in individual shells around imol
          m_frac_comp : fractional composition of lipid types around imol
    """

    # generate shells
    member_shell: list[set[int]] = []

    # shell index
    id_shell = [-1 for _ in range(nmol)]

    # set for all molecules
    all_list: set[int] = set()
    for i in range(0, nmol):
        all_list.add(i)

    # set for already assigned molecules
    assigned_list = set()

    # iteratively generate shell member
    # 0th shell: self
    ishell = 0
    member_shell.append(set())
    member_shell[ishell].add(imol)
    id_shell[imol] = ishell
    assigned_list.add(imol)

    # generate member_shell[ishell+1] from member_shell[ishell]
    # This won't work properly at large shells: So limit the max_shell
    #
    # 1. Set of contactMol for member_shell[ishell]
    #      outer_shell + current_shell (+ inner_shell (for ishell >0))
    #
    # 2. subtract current_shell members (& inner_shell members)
    for ishell in range(0, max_shell):
        current_shell_memb: set[int] = member_shell[ishell]
        tset: set[int] = set()
        for elem in current_shell_memb:
            # add contactMol[elem] to the next shell members
            tset = tset.union(contactMol[elem])
        # trim tset by subtracting current shell members
        tset = tset - current_shell_memb
        # inner shell members to be subtracted for ishell >0
        if ishell > 0:
            tset = tset - member_shell[ishell-1]
        # outer shell members
        member_shell.append(tset)  # member_shell[ishell+1]

    # setup number and fractional compositions for individual shells
    m_numb_comp = np.zeros([max_shell+1, ntype], dtype=float)
    m_frac_comp = np.zeros([max_shell+1, ntype], dtype=float)

    # calculate number and fractional compositions up to max_shell-th shell
    for i in range(0, max_shell+1):
        tmemb_shell = member_shell[i]
        for elem in tmemb_shell:
            j = id_type[elem]       # mol type index
            m_numb_comp[i][j] += 1.0  # update counter
        norm = np.sum(m_numb_comp[i])
        if norm > 0.0:
            m_frac_comp[i] = m_numb_comp[i]/norm

    return m_numb_comp, m_frac_comp


def update_shell_comp(ntype: int, max_shell: int,
                      t_numb_comp: NDFloat64, t_frac_comp: NDFloat64,
                      tnumb_comp: NDFloat64, tfrac_comp: NDFloat64) \
        -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Update raw shell compositions
    ----------
    NOTE:
        update shell-wise compositions from data for a molecule
        returns updated max_shell & accumulated (i.e., un-normalized) shell compositions

    input
          ntype      : number of unique molecule types
          max_shell  : maximum shell index (set to be the same for analysis)
          t_numb_comp: t_numb_comp[frame (,leaflet), type, shell, type]
          t_frac_comp: t_frac_comp[frame (,leaflet), type, shell, type]
          tnumb_comp : shell-wise number compositions of lipid types around a molecule
          fnumb_comp : shell-wise fractional compositions of lipid types around a molecule

    ouput
          t_numb_comp: updated t_numb_comp[frame (,leaflet), type, shell, type]
          t_frac_comp: updated f_numb_comp[frame (,leaflet), type, shell, type]
    """

    # update shell compositions from leaflets at the frame, iframe
    for ishell in range(0, max_shell+1):
        for k in range(0, ntype):
            t_numb_comp[ishell][k] += tnumb_comp[ishell][k]
            t_frac_comp[ishell][k] += tfrac_comp[ishell][k]

    return t_numb_comp, t_frac_comp


def normalize_raw_comp(array1: NDFloat64, array2: NDFloat64,
                       array3: Sequence[float] | NDFloat64) \
        -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Normalize shell-wise compositions
    ----------

    input
          array1: numb_comp[frame (,leaflet), type1, shell, type2]
          array2: frac_comp[frame (,leaflet), typd1, shell, type2]
          array3: number of type1 (normalization factor); dim. = [frame (,leaflet), type1]

    output
          array1: normalized numb_comp
          array2: normalized frac_comp
    """

    ntype = len(array1)     # len(array2)     # len(array3)
    # len(array2[0])  # - number of shells to get normalized compositions
    # nshell = len(array1[0]) - should be different between types

    for i in range(0, ntype):
        tnorm = array3[i]
        if tnorm > 0.0:
            nshell = len(array1[i])  # shell number for i-th type
            for j in range(0, nshell):
                for k in range(0, ntype):
                    array1[i][j][k] /= tnorm
                    array2[i][j][k] /= tnorm

    return array1, array2


# TBD
# Clustering analysis using Voronoi tessellation
# contact info. from the above functions will be used
# It is not clear --- LEAVE as ongoing but do not include yet
# -------
# #
# combine N.N. VCs infor sorted by weightVCs
# #
# ----------
# input
# contactVC    : list of NNVCs for ind. VC
# borderVC: list array of border lengths btw. ind. VC and its NN VCs.
# distVC  : list array of distances btw. ind. VC and its NN VCs
# weightVC: list array of weights of NNVCs for ind. VC
# output
# NNVCinfo: list of array of NN data for individual VCs
# #
# def combine_NNVCinfo(contactVC,distVC,borderVC,weightVC):
# define a data type
# dtype = [ ('ivc',int) , ('weight',float) , ('border',float) , ('dist',float)]
##
# sorted NN info; structured list array !
# NNVCinfo = []
##
# loop over VCs
# for i in range(0,nprim):
# nnvc = contactVC[i]
# weight = weightVC[i]
# border = borderVC[i]
# dist   = distVC[i]
# num_nnvc = len(nnvc)
##
# tarray = []
# for j in range(0,num_nnvc):
# tivc = nnvc[j]
# twght= weight[j]
# tbrd = border[j]
# tdist= dist[j]
# tarray.append([tivc,twght,tbrd,tdist])
##
# tnparray = np.array(tarray, dtype=dtype)
##
# sort tnparray along 'weight' - ascending order
# tnnvcinfo = np.sort(tnparray, order='weight',stable=True)
# NNVCinfo.append( tnnvcinfo[::-1])
##
# return NNVCinfo


# --------------------
# Functions for area analysis
# --------------------

def calc_vr_area(ivc: int, vorn: 'Voronoi') -> float:
    """
    ----------
    Calculate area of a Voronoi region
    ----------

    input
          ivc: voronoi center index
          ivt: vertices for the region (polygon) surrounding the VC
    output
          area: area of the voronoi region for the VC,iv
    """

    ipr = vorn.point_region[ivc]  # index of region for the VC, ivc
    ivt = vorn.regions[ipr]      # indices of vertices for the voronoi region

    area = 0.0
    nvert = len(ivt)  # number of vertices
    if nvert > 0:
        tinx = ivt[0]  # staring index of vertices
        vt1 = vorn.vertices[tinx]
        vt0 = vt1
        # loop over vertices
        for m in range(1, nvert):
            tinx = ivt[m]
            vt2 = vorn.vertices[tinx]
            area += np.cross(vt1, vt2)
            vt1 = vt2  # update vertex coordinates
        # last pair
        area += np.cross(vt1, vt0)
        area = math.fabs(0.5*area)

    return area


def calc_mol_area(mol_ivc: list['Voronoi'], vorn: 'Voronoi') -> float:
    """
    ----------
    Calculate molecular area
    ----------

    input
          mol_ivc: VC list for a chosen molecule
          vorn   : output vorn=Voronoi(vcent); Voronoi() from scipy

    output
          area: area of the molecule
    """

    area = 0.0
    for ivc in mol_ivc:
        area += calc_vr_area(ivc, vorn)
    return area


def calc_mol_area_all(nmol: int, mol_ivc: list['Voronoi'],
                      vorn: 'Voronoi') -> NDFloat64:
    """
    ----------
    Calculate area of all molecules
    ----------

    input
          nmol   : number of molecules
          mol_ivc: stariting VC index for individual molecules
          vorn   : vorn=Voronoi(vcent)

    output
          area   : array of areas of all molecules
    """

    # set area array
    area = np.zeros([nmol], dtype=float)

    for i in range(0, nmol):
        area[i] = calc_mol_area(mol_ivc[i], vorn)

    return area


def calc_apl(area: NDFloat64, id_type: list[int], ntype: int) -> tuple[NDFloat64, NDFloat64]:
    """
    ----------
    Calculate area per lipids for individual molecule types
    ----------

    input
          area   : area of individual molecules
          id_type: molecule type index of individual molecules
          ntype  : number of unique molecule types

    output
          apl    : component area per lipid
          ncnt   : number of lipids for each types
    """

    # component APL
    apl = np.zeros([ntype], dtype=float)

    # number of each lipid type
    ncnt = np.zeros([ntype], dtype=float)

    # number of molecules
    nmol = len(area)

    for i in range(0, nmol):
        # id_type_vc[mol_ivc[i][0]]      # lipid-type indices
        itype = id_type[i]
        # print(f'# id_type = {itype}')
        apl[itype] += area[i]
        ncnt[itype] += 1.0

    # get component APL
    for i in range(0, ntype):
        tsum = ncnt[i]
        if tsum > 0.0:
            apl[i] /= tsum

    return apl, ncnt


def calculate_ave_and_std(nside: int, ntype: int, apl: NDFloat64, ncomp: NDInt64,
                          outt: OutputFileType = 'outl') -> tuple[NDFloat64, NDFloat64, NDFloat64]:
    """
    ----------
    Calculate average and std of APLs
    ----------

    input
          nside : number of leaflets
          ntype : number of qunique lipid types
          apl   : time series of component APL
          ncomp : counts for the component APL
          outt  : output type. outl (leaflets) or outb (bilayer)

    output
          aapl  : average component APL
          sapl  : std of component APL
          anum  : total counts for the component APL
    """
    def do_small() -> tuple[NDFloat64, NDFloat64, NDFloat64]:
        aapl = np.zeros([ntype], dtype=float)  # mean
        sapl = np.zeros([ntype], dtype=float)  # std
        anum = np.zeros([ntype], dtype=float)  # number
        for j in range(0, ntype):
            tot_weight = np.sum(ncomp[:, j], axis=0)  # total weight
            anum[j] = tot_weight
            if tot_weight > 0.0:
                weight = ncomp[:, j]/tot_weight
                tapl = apl[:, j]  # apl time series

                ave = np.average(tapl, weights=weight)
                std = np.average((tapl-ave)**2, weights=weight)
                std = math.sqrt(std)

                aapl[j] = ave
                sapl[j] = std

        return aapl, sapl, anum

    def do_big() -> tuple[NDFloat64, NDFloat64, NDFloat64]:
        aapl = np.zeros([nside, ntype], dtype=float)  # mean
        sapl = np.zeros([nside, ntype], dtype=float)  # std
        anum = np.zeros([nside, ntype], dtype=float)  # number
        for i in range(0, nside):
            for j in range(0, ntype):
                tot_weight = np.sum(ncomp[:, i, j], axis=0)  # total weight
                anum[i, j] = tot_weight
                if tot_weight > 0.0:
                    weight = ncomp[:, i, j]/tot_weight
                    tapl = apl[:, i, j]  # apl time series

                    ave = np.average(tapl, weights=weight)
                    std = np.average((tapl-ave)**2, weights=weight)
                    std = math.sqrt(std)

                    aapl[i, j] = ave
                    sapl[i, j] = std

        return aapl, sapl, anum

    # APL:  WEIGHTED MEAN AND STD
    if outt == 'outb':
        return do_small()
    return do_big()


# --------------------
# Other functions
# --------------------

def calc_com_all(nmol: int, mol_ivc: List2D[int], vorn: 'Voronoi', box: NDFloat64) -> NDFloat64:
    """
    ----------
    Calculate XY-COM of molecules from associated VCs
    ----------
    NOTE: periodicity is considered.

    input
          nmol   : number of molecules
          mol_ivc: stariting VC index for individual molecules
          vorn   : vorn=Voronoi(vcent)

    output
          com    : array of XY-COM of all molecules
    """

    # set com array
    com = np.zeros([nmol, 2], dtype=float)

    hbox = box/2.0

    # loop over molecules
    for i in range(0, nmol):
        nsite = len(mol_ivc[i])
        norm = float(nsite)

        p0 = vorn.points[mol_ivc[i][0]]  # first point
        tcom = 0 * p0
        if nsite > 1:
            for j in range(1, nsite):
                p1 = vorn.points[mol_ivc[i][j]]
                disp = p1 - p0  # displacement
                disp = disp - np.sign(np.trunc(disp/hbox))*box
                tcom += disp
        tcom /= norm
        tcom += p0

        # make tcom be in the box
        disp = tcom - box/2
        disp = disp - np.sign(np.trunc(disp/hbox))*box
        tcom = disp + box/2

        com[i] = tcom

    return com
