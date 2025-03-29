#!/usr/bin/python
import sys,math
import numpy as np 

## #--- The following are hard set for voronoi tessellation in 2D
## dim=2              # 2 dimensions
## nimage=int(3**dim) # number of total primary+images for 2D

def connect_mol_and_vc(u,nmol,nmol_type,ag):
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
  for i in range(0,nmol): npt += len(ag[i])
  print(f'# total number of VCs = {npt}')

  id_mol_vc  = np.zeros([npt],dtype=int)
  mol_ivc    = []

  iatom = 0
  for i in range(0,nmol):
    mol_ivc.append([])
    natom = len(ag[i])
    for j in range(0,natom):
      mol_ivc[i].append(iatom)
      id_mol_vc[iatom] = i
      iatom += 1
  
  return (mol_ivc,id_mol_vc)



def image_translation(box): 
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
  xtla,xtlb =box[0],box[1]

  # translation sequence (primary, ..., last image)
  # Hard set operator for image translation in 2D
  operator = np.array([[0,0],  [1,0],  [1,1],\
                       [0,1],  [-1,1],[-1,0],\
                       [-1,-1],[0,-1], [1,-1]], dtype=float)
  noperator = len(operator) # ! should be equal to nimage

  # translation vectors
  tran=np.zeros([noperator,2],dtype=float)
  for i in range(0,noperator):
    np.copyto(tran[i,:],operator[i,:]*box)
  
  return (tran)



def prepare_voronoi_centers(ag,box,dim,nimage): 
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
  nmol = len(ag) # atom groups for individual molecules 
  for i in range(0,nmol): nprim += len(ag[i])
  nvorn=int(nimage*nprim)    # total VCs in 2D (primary + images)

  sout =f'# DIMENSION= {dim} NIMAGE= {nimage}\n'
  sout+=f'# VORONOI CENTER INFO.\n'
  sout+=f'       nprim= {nprim} nvorn= {nvorn}\n'
  print(sout)

  # initialize VC array
  vcent=np.zeros([nimage*nprim,dim],dtype=float)
  tvcen=np.zeros([nprim,dim],dtype=float)

  iatom = 0
  for i in range(0,nmol):
    natom = len(ag[i]) # number of atoms
    ist,ied = (iatom,iatom+natom) # start/end indices in tvcen
    tpos = ag[i].positions
    np.copyto(tvcen[ist:ied],tpos[:,0:dim]) # copy xy coordinates to tvcen
    iatom += natom

  # translation vector for all images
  tran = image_translation(box)

  # loop over images
  for iimg in range(0,nimage):

    tmpvc=tvcen + tran[iimg] # translation

    # start & stop indices in vcent array
    ist,ied = int(iimg*nprim),int((iimg+1)*nprim)

    # copy image coordinates to vcent array
    np.copyto(vcent[ist:ied],tmpvc)
    # print(f'# iimg={iimg}, xtla={box[0]} xtlb={box[1]} vcent[ist]= {vcent[ist]}')

  return (nprim,nvorn,vcent)



#--------------------
# Functions for CONTACT analysis
#--------------------

def generate_vc_contact_info(vorn,nprim,nimage,box):
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
  contactVC, borderVC, distVC, weightVC = [],[],[],[]
  for i in range(0,nprim):
    contactVC.append([])
    borderVC.append([])
    distVC.append([])
    weightVC.append([])

  # search for vorn.ridge_points
  nridges=len(vorn.ridge_points)
  for i in range(0,nridges):
    inx1,inx2 = vorn.ridge_points[i,0], vorn.ridge_points[i,1]

    # check validity
    if (inx1 < 0 or inx2 < 0):
      print(f'# construct N.N. VC: Error - index range: {inx1} {inx2}')
      sys.exit(1)
    # generate only for the VCs in the primary image
    if (inx1 >= nprim and inx2 >= nprim): continue

    # inx1 < nprim or inx2 < nprim
    # calc border length btw. VCs
    vinx = vorn.ridge_vertices[i] # indices of the voronoi vertices for the ridge
    vert1 = vorn.vertices[vinx[0]]
    vert2 = vorn.vertices[vinx[1]]
    tborder = np.linalg.norm(vert1 - vert2)
    # calc dist btw. VCs
    tdist = np.linalg.norm(vorn.points[inx1]-vorn.points[inx2])

    # 1.  add indices to contactVC & update distVC and borderVC
    if (inx1 < nprim):
      contactVC[inx1].append(inx2)
      borderVC[inx1].append(tborder)
      distVC[inx1].append(tdist)
      weightVC[inx1].append(tborder/tdist)

    if (inx2 < nprim):
      contactVC[inx2].append(inx1)
      borderVC[inx2].append(tborder)
      distVC[inx2].append(tdist)
      weightVC[inx2].append(tborder/tdist)

  # convert into numpy array
  for i in range(0,nprim):
    tarray = contactVC[i];     contactVC[i]     = np.array(tarray);
    tarray = borderVC[i]; borderVC[i] = np.array(tarray);
    tarray = distVC[i];   distVC[i]   = np.array(tarray);
    tarray = weightVC[i]; weightVC[i] = np.array(tarray);

  return (contactVC,borderVC,distVC,weightVC)



def generate_mol_contact_info(mol_ivc,contactVC,id_mol_vc,nprim):
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
  nmol  = len(mol_ivc) 

  # empyt list output array
  contactMol = []

  # Loop over molecules
  for i in range(0,nmol):
    # collect all contactVCs for molecule, i
    tset = set()  # empty set for the contactVCs for molecule, i
    nsite=len(mol_ivc[i]) 
    for elem1 in mol_ivc[i]:
      for elem2 in contactVC[elem1]:
        tid = id_mol_vc[elem2 % nprim] # molecule index 
        if (tid != i): tset.add(tid)

    # append to the output array
    contactMol.append(tset)
   
  return (contactMol)



def calculate_contact(nmol,ntype,nmol_type,mol_ivc,id_type,contactMol):
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
  ncontact = np.zeros([ntype,ntype],dtype=float)
  fcontact = np.zeros([ntype,ntype],dtype=float)

  # loop over mols
  for i in range(0,nmol):
    inx1 = id_type[i] 
    tnn_mol = contactMol[i]
    for j in tnn_mol:
      inx2 = id_type[j] 
      ncontact[inx1][inx2] += 1.0

  # division by the number of components
  for i in range(0,ntype):
    norm = nmol_type[i]
    if (norm > 0.0):
      ncontact[i,:]/= norm

  # calculate fcontact
  for i in range(0,ntype):
    tnorm = np.sum(ncontact[i,:])
    if (tnorm > 0.0):
      fcontact[i] = ncontact[i]/tnorm

  return (ncontact,fcontact)

#--------------
#
# Update contact for bilayer from leaflet data
#
#--------------
def update_contact_bilayer(framenum,ntype,ncontact,fcontact,tnmol_type):
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

  for i in range(0,framenum):
    # update ncontact
    for j in range(0,ntype):
      tnorm = tnmol_type[i,j] # it's component should not be zero.
      ncontact[i,j,:] /= tnorm
    # update fcontact
    for j in range(0,ntype):
      tnorm = np.sum(ncontact[i,j,:])
      if (tnorm > 0.0): fcontact[i,j,:] = ncontact[i,j,:]/tnorm

  return ncontact,fcontact



#--------------------
# Functions for Shell analysis
#--------------------

def setup_shell_comp_arrays(framenum,nside,ntype,qb):
  """
  ----------
  Setup raw shell-wise composition arrays
  ----------

  input
        framenum : number of frames for the analysis
        nside    : 2 (bilayer)
        ntype    : number of unique molecule types
        qb       : Options for analysis; True (bilayer)/False (leaflet)
  
  output
        max_shell: maximum number of shells around a given molecule type
        t_numb_comp: list array of time series of shell-wise number compositions 
        t_frac_comp: list array of time series of shell-wise fractional compositions
  """

  max_shell = []
  t_numb_comp,t_frac_comp = [],[]
  if (qb is True):
    for i in range(0,framenum):         
      max_shell.append([])
      t_numb_comp.append([])
      t_frac_comp.append([])
      for j in range(0,ntype):
        max_shell[i].append(100) # give a large number & update after analysis
        t_numb_comp[i].append([])
        t_frac_comp[i].append([])
  else:   
  # if (qb is False):
    for i in range(0,framenum):
      max_shell.append([])
      t_numb_comp.append([])
      t_frac_comp.append([])
      for j in range(0,nside):
        max_shell[i].append([])
        t_numb_comp[i].append([])
        t_frac_comp[i].append([])
        for k in range(0,ntype):
          max_shell[i][j].append(100) # give a large number & update after analysis
          t_numb_comp[i][j].append([])
          t_frac_comp[i][j].append([])
  return (max_shell,t_numb_comp,t_frac_comp)



def calculate_shell_comp(imol,nmol,contactMol,ntype,mol_ivc,id_type):
  """
  ----------
  Calculate shell-wise number and fractional compositions around a molecule
  ----------

  input
        imol      : molecule index to generate shell arount it
        nmol      : total number of molecules
        contactMol: contacting molecule index list for individual molecule
        mol_ivc   : voronoi cell list for individual molecules
        id_type   : lipid type index for individual molecule
  
  output
        max_shell   : max shell index
        numb_comp   : numbers of lipid types in individual shells
        frac_comp   : fractional composition of lipid types
  """

  # generate shells
  member_shell = []    

  # shell index
  id_shell = [] ; [id_shell.append(-1) for i in range(0,nmol)];

  # set for already assigned molecules
  assigned_list = set()
  
  # iteratively generate shell member
  # 0th shell: self
  ishell = 0
  member_shell.append(set())
  member_shell[ishell].add(imol)
  id_shell[imol] = ishell
  assigned_list.add(imol)

  # read contactMol from the shell & subtract already assigned molecule indices
  while (True):
    tset = set() # temp set for the next shell
    for elem1 in member_shell[ishell]:
      tset=tset.union(contactMol[elem1]) # all contactMol for the shell member
    # print(f'tset',tset)
    # subtract the previously assigned molecules
    tset = tset - assigned_list
    if (len(tset) == 0):
      max_shell = ishell 
      # print(f'# max_shell = {max_shell}');
      break;
    else:
      # update the next shell members
      ishell += 1
      member_shell.append(set().union(tset))
      for elem2 in member_shell[ishell]: id_shell[elem2] = ishell
      assigned_list=assigned_list.union(tset)
      ## print(f'shell {ishell}',member_shell[ishell])
      ## print(f'assigned_list', assigned_list)

  # setup number and fractional compositions for individual shells
  numb_comp = []
  frac_comp = []
  for i in range(0,max_shell+1):
    numb_comp.append([])
    frac_comp.append([])
    numb_comp[i] = np.zeros([ntype],dtype=float)
    frac_comp[i] = np.zeros([ntype],dtype=float)
    # for j in range(0,ntype):
    #   frac_comp[i].append(0.0) 

  # calculate number and fractional compositions
  for i in range(0,max_shell+1):
    tmemb_shell = member_shell[i]
    for elem in tmemb_shell:
      j = id_type[elem]       # mol type index
      numb_comp[i][j] += 1.0  # update counter
    norm = np.sum(numb_comp[i])
    if (norm > 0.0): 
      frac_comp[i] = numb_comp[i]/norm

  return (max_shell,numb_comp,frac_comp)



def update_shell_comp(ntype,max_shell,t_numb_comp,t_frac_comp,\
                      tmax_shell,tnumb_comp,tfrac_comp):
  """
  ----------
  Updaet raw shell compositions
  ----------
  NOTE:
      update shell-wise compositions from data for a molecule
      returns updated max_shell & accumulated (i.e., un-normalized) shell compositions

  input
        ntype      : number of unique molecule types
        max_shell  : max_shell[frame (,leaflet), type]
        t_numb_comp: t_numb_comp[frame (,leaflet), type]
        t_frac_comp: t_frac_comp[frame (,leaflet), type]
        tmax_shell : max_shell for a given molecule
        tnumb_comp : shell-wise number compositions of lipid types around a molecule
        fnumb_comp : shell-wise fractional compositions of lipid types around a molecule
  
  ouput
        max_shell  : updated max_shell[frame (,leaflet), type]
        t_numb_comp: updated t_numb_comp[frame (,leaflet), type, shell, type]
        t_frac_comp: updated f_numb_comp[frame (,leaflet), type, shell, type]
  """

  # - update shell compositions from leaflets at the frame, iframe
  # - expand partially setup composition arrays for the frame
  # - then update compositions (by adding frame data)

  # update time series data
  tmax_shell0 = max_shell
  # min of max_shell 
  if (tmax_shell < tmax_shell0):
    max_shell = tmax_shell
  # expand t_num_comp/t_frac_comp to next level [shell]

  # Check size of array at shell level
  # and expand upto tmax_shell+1  
  # then initialize shell compositions 
  len1 = len(t_numb_comp)
  len2 = len(t_frac_comp)
  if (len1 != len2):
    print(f'# incosistent array size between numb & frac compositions!')
    sys.exit(1)
  if (len1 < tmax_shell+1):
    for ishell in range(len1,tmax_shell+1):
      t_numb_comp.append([])
      t_frac_comp.append([])
      for k in range(0,ntype):
        t_numb_comp[ishell].append(0.0)
        t_frac_comp[ishell].append(0.0)

  # update compositions
  for ishell in range(0,tmax_shell+1):
    for k in range(0,ntype):
      t_numb_comp[ishell][k] += tnumb_comp[ishell][k]
      t_frac_comp[ishell][k] += tfrac_comp[ishell][k]

  return (max_shell,t_numb_comp,t_frac_comp)



def normalize_raw_comp(array1,array2,array3):
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
  nshell = len(array1[0]) # len(array2[0])  # - number of shells to get normalized compositions
  
  for i in range(0,ntype):
    tnorm = array3[i]
    if (tnorm > 0.0):
      for j in range(0,nshell):
        for k in range(0,ntype):
          array1[i][j][k] /= tnorm
          array2[i][j][k] /= tnorm
  
  return (array1,array2)



def set_nshell(max_shell,framenum,nside,ntype,qb):
  """
  ----------
  Set maximum shell number for statistics and outputs
  ----------

  input
        max_shell: array of maximum number of shells; dim = [frame (,leaflet), type]
        framenum : number of frames for the analysis
        nside    : number of leaflets, 2
        ntype    : number of unique molecule types
        qb       : Option for the analysis; True (bilayer)/ False (leaflet)
  
  output
        nshell   : maximum number of shells for statistics and outputs
  """

  nshell = 100 # global value
  if (qb is True):
    for i in range(0,framenum):
      for j in range(0,ntype):
        tmax_shell = max_shell[i][j]
        # print(tmax_shell)
        if (nshell > tmax_shell): nshell = tmax_shell
  else:
  # if (qb is False):
    for j in range(0,nside):
      for i in range(0,framenum):
        for k in range(0,ntype):
          tmax_shell = max_shell[i][j][k]
          # print(tmax_shell)
          if (nshell > tmax_shell): nshell = tmax_shell

  print(f'# nshell= {nshell}') # min. max shell index for all lipid types.
  return (nshell)



def convert_to_ndarray(t_numb_comp,t_frac_comp,framenum,nside,ntype,nshell,qb):
  """
  ----------
  Convert list arrays of time series to ndarrays for simpler statistics
  ----------

  input
        t_numb_comp: time series of shell-wise number compositions
        t_frac_comp: time series of shell-wise fractional compositions
        framenum   : number of frames for the analysis
        nside      : number of leaflets, 2
        ntype      : number of unique molecule types
        nshell     : maximum shell number for statistics and outputs
        qb         : Option for the analysis; True (bilayer)/ False (leaflet)
  
  output
        numb_comp  : ndarray of time series of shell-wise number compositions
        frac_comp  : ndarray of time series of shell-wise fractional compositions
  """

  if (qb is True):
    # Statistics upto min_max_shell
    numb_comp = np.zeros([framenum,ntype,nshell+1,ntype],dtype=float)
    frac_comp = np.zeros([framenum,ntype,nshell+1,ntype],dtype=float)
    # sum up value & get mean and ave
    for i in range(0,framenum):
      for j in range(0,ntype):
        for k in range(0,nshell+1):
          for m in range(0,ntype):
            numb_comp[i,j,k,m] = t_numb_comp[i][j][k][m]
            frac_comp[i,j,k,m] = t_frac_comp[i][j][k][m]
  else:  # if (qb is False):
    # Statistics upto min_max_shell
    numb_comp = np.zeros([framenum,nside,ntype,nshell+1,ntype],dtype=float)
    frac_comp = np.zeros([framenum,nside,ntype,nshell+1,ntype],dtype=float)
    # sum up value & get mean and ave
    for i in range(0,framenum):
      for j in range(0,nside):
        for k in range(0,ntype):
          for m in range(0,nshell+1):
            for n in range(0,ntype):
              numb_comp[i,j,k,m,n] = t_numb_comp[i][j][k][m][n]
              frac_comp[i,j,k,m,n] = t_frac_comp[i][j][k][m][n]

  return (numb_comp,frac_comp)

## TBD
## Clustering analysis using Voronoi tessellation
##	contact info. from the above functions will be used
## #### It is not clear --- LEAVE as ongoing but do not include yet
## #-------
## #
## # combine N.N. VCs infor sorted by weightVCs
## #
## #----------
## # input
## #	contactVC    : list of NNVCs for ind. VC
## #	borderVC: list array of border lengths btw. ind. VC and its NN VCs. 
## #	distVC  : list array of distances btw. ind. VC and its NN VCs
## #	weightVC: list array of weights of NNVCs for ind. VC
## # output
## #	NNVCinfo: list of array of NN data for individual VCs
## #
## def combine_NNVCinfo(contactVC,distVC,borderVC,weightVC):
##   # define a data type
##   dtype = [ ('ivc',int) , ('weight',float) , ('border',float) , ('dist',float)]
## 
##   # sorted NN info; structured list array !
##   NNVCinfo = []
## 
##   # loop over VCs
##   for i in range(0,nprim):
##     nnvc = contactVC[i]
##     weight = weightVC[i]
##     border = borderVC[i]
##     dist   = distVC[i]
##     num_nnvc = len(nnvc)
## 
##     tarray = []
##     for j in range(0,num_nnvc):
##       tivc = nnvc[j]
##       twght= weight[j]
##       tbrd = border[j]
##       tdist= dist[j]
##       tarray.append([tivc,twght,tbrd,tdist]) 
## 
##     tnparray = np.array(tarray, dtype=dtype)
## 
##     # sort tnparray along 'weight' - ascending order
##     tnnvcinfo = np.sort(tnparray, order='weight',stable=True)
##     NNVCinfo.append( tnnvcinfo[::-1])
## 
##   return (NNVCinfo)



#--------------------
# Functions for area analysis
#--------------------

def calc_vr_area(ivc,vorn):
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

  ipr=vorn.point_region[ivc] # index of region for the VC, ivc
  ivt=vorn.regions[ipr]      # indices of vertices for the voronoi region

  area=0.0
  nvert=len(ivt) # number of vertices
  if(nvert > 0): 
    tinx=ivt[0]  # staring index of vertices
    vt1=vorn.vertices[tinx] ; vt0=vt1;
    # loop over vertices
    for m in range(1,nvert): 
      tinx=ivt[m]
      vt2=vorn.vertices[tinx]
      area += np.cross(vt1,vt2)
      vt1 = vt2 # update vertex coordinates
    # last pair
    area += np.cross(vt1,vt0)
    area = math.fabs(0.5*area)
  
  return area



def calc_mol_area(mol_ivc,vorn):
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

  area=0.0;
  for ivc in mol_ivc:
    area += calc_vr_area(ivc,vorn)
  return area



def calc_mol_area_all(nmol,mol_ivc,vorn):
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
  area=np.zeros([nmol],dtype=float)

  for i in range(0,nmol):
    area[i]=calc_mol_area(mol_ivc[i],vorn)

  return area



def calc_apl(area,id_type,ntype):
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
  apl = np.zeros([ntype],dtype=float)

  # number of each lipid type 
  ncnt = np.zeros([ntype],dtype=float)

  # number of molecules
  nmol=len(area)

  for i in range(0,nmol):
    itype=id_type[i] # id_type_vc[mol_ivc[i][0]]      # lipid-type indices 
    # print(f'# id_type = {itype}')
    apl[itype] += area[i]
    ncnt[itype] += 1.0

  # get component APL
  for i in range(0,ntype):
   tsum=ncnt[i]
   if (tsum > 0.0): apl[i] /= tsum

  return (apl,ncnt)



def calculate_ave_and_std(nside,ntype,apl,ncomp,qb):
  """
  ----------
  Calculate average and std of APLs
  ----------

  input
        nside : number of leaflets
        ntype : number of qunique lipid types
        apl   : time series of component APL
        ncomp : counts for the component APL
        qb    : True (bilayer)/ False (leaflet)
  
  output
        aapl  : average component APL
        sapl  : std of component APL
        anum  : total counts for the component APL
  """

  # APL:  WEIGHTED MEAN AND STD
  if (qb is True):
    aapl = np.zeros([ntype],dtype=float) # mean
    sapl = np.zeros([ntype],dtype=float) # std
    anum = np.zeros([ntype],dtype=float) # number
    for j in range(0,ntype):
      tot_weight = np.sum(ncomp[:,j],axis=0)   #  total weight
      anum[j] = tot_weight
      if (tot_weight > 0.0):
        weight = ncomp[:,j]/tot_weight
        tapl = apl[:,j]                        #  apl time series

        ave=np.average(tapl,weights=weight)
        std=np.average((tapl-ave)**2, weights=weight)
        std=math.sqrt(std)

        aapl[j] = ave
        sapl[j] = std
  else: # if (qb is False):
    aapl = np.zeros([nside,ntype],dtype=float) # mean
    sapl = np.zeros([nside,ntype],dtype=float) # std
    anum = np.zeros([nside,ntype],dtype=float) # number
    for i in range(0,nside):
      for j in range(0,ntype):
        tot_weight = np.sum(ncomp[:,i,j],axis=0)   #  total weight
        anum[i,j] = tot_weight
        if (tot_weight > 0.0):
          weight = ncomp[:,i,j]/tot_weight
          tapl = apl[:,i,j]                        #  apl time series

          ave=np.average(tapl,weights=weight)
          std=np.average((tapl-ave)**2, weights=weight)
          std=math.sqrt(std)

          aapl[i,j] = ave
          sapl[i,j] = std

  return (aapl,sapl,anum)



#--------------------
# Other functions
#--------------------

def calc_com_all(nmol,mol_ivc,vorn,box):
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
  com=np.zeros([nmol,2],dtype=float)

  hbox = box/2.0

  # loop over molecules
  for i in range(0,nmol):
    nsite = len(mol_ivc[i])
    norm = float(nsite)

    p0 = vorn.points[mol_ivc[i][0]] # first point
    tcom = 0 * p0
    if (nsite > 1):
      for j in range(1,nsite):
        p1 = vorn.points[mol_ivc[i][j]]
        disp = p1 - p0 # displacement
        disp = disp - np.sign(np.trunc(disp/hbox))*box
        tcom += disp
    tcom /= norm
    tcom += p0

    # make tcom be in the box
    disp = tcom - box/2
    disp = disp - np.sign(np.trunc(disp/hbox))*box
    tcom = disp + box/2

    com[i]=tcom

  return (com)

