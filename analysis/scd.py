#!/usr/bin/python
import argparse
from typing import Optional
import os,sys,re,math
import stanalyzer.bin.stanalyzer as sta
import numpy as np       
import MDAnalysis as mda 
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup
from . import leaflet_util as myleaflet
from . import mol_atomgroup_util as mymol
from . import scd_util as myscd

ANALYSIS_NAME = 'scd'


#--- The following are hard set for membrane analysis
nside=2     # up/dn
sside=["up","dn"]



def process_args(sel,split):
  """
  ----------
  Process arguments
  ----------
  """

  selection = re.split(';|,', f'{sel:s}') # selections for individual lipid types
  ntype = len(selection)                  # number of unique lipid types
  for i in range(0,ntype): selection[i]=selection[i].strip()

  # Process selection strings to extract information
  sel_type  = [] # selection type
  name_type = [] # name of lipid type
  ic_name   = [] # starting carbon names for ind. tail in ind. lipid type
  nchain    = [] # number of carbon chains in ind. lipid types
  ref_name  = [] # reference atom for leaflet assignment
  for i in range(0,ntype):
    tmps=selection[i].split() 
    sel_type.append(tmps[0])   # segid/resname/moleculetype
    name_type.append(tmps[1])  # PROA/PRO*/DSPC/...
    ic_name.append([])         # 
    ref_name.append([])

    # process remaining strings
    for j in range(2,len(tmps)):
      tmps_check = tmps[j].strip("(").strip(")")
      tstr = tmps_check.lower()
      if(len(tstr) == 0): continue
      if(tstr == "name" or tstr == "and" or tstr == "or" or tstr == "not"):
        continue
      if (j < len(tmps) - 1):
        ic_name[i].append(tmps_check) # append starting carbon name 
      else:
        ref_name[i].append(tmps_check) # last atom; ref atom
    # update nchain
    nchain.append(len(ic_name[i]))

  split  = re.split(';|,',f'{split:s}')
  nsplit = len(split) 
  for i in range(0,nsplit): split[i]=split[i].strip()
  qsplit = []
  if (nsplit < ntype): # add more qsplit options
    for i in range(0,nsplit):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)
    for i in range(nsplit,ntype): qsplit.append(True) # default values
  else: #  get args.split upto ntype
    qsplit = []
    for i in range(0,ntype):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)
  return (selection,ntype,sel_type,name_type,ref_name,ic_name,nchain,qsplit)



def write_ave_std_leaflet(nside,ntype,name_type,sside,nchain,carbons,ncarbons,array1,array2,odir):
  # sside : leafelt names
  # array1: average
  # array2: std
  # odir  : output path
  for i in range(0,nside):
    side=sside[i]
    # print for lipid types
    for j in range(0,ntype):
      ntail = nchain[j]
      for k in range(0,ntail):
       sout = f'# {side}: {name_type[j]}  chain{k} SCD\n'
       # loop over carbon 
       nc = ncarbons[j][k]
       for m in range(0,nc):
         sout += f' {carbons[j][k][m]:6s} {array1[i][j][k][m]:10.5f} {aarray2[i][j][k][m]:10.5f}\n'
       print(sout)
       fout=f'{odir}/{side}_{name_type[j].lower()}_chain{k}.plo'
       f = open (fout,'w'); f.write(sout) ; f.close();
  return



def write_ave_std_bilayer(ntype,name_type,nchain,carbons,ncarbons,array1,array2,odir):
  # array1: average
  # array2: std
  # odir  : output path
  side="bilayer"
  # print for lipid types
  for j in range(0,ntype):
    ntail = nchain[j]
    for k in range(0,ntail):
     sout = f'# {side}: {name_type[j]} chain{k} SCD\n' 
     # loop over carbon 
     nc = ncarbons[j][k]
     for m in range(0,nc):
       sout += f' {carbons[j][k][m]:6s} {array1[j][k][m]:10.5f} {array2[j][k][m]:10.5f}\n'
     print(sout)
     fout=f'{odir}/{side}_{name_type[j].lower()}_chain{k}.plo'
     f = open (fout,'w'); f.write(sout) ; f.close();
  return



def write_time_series_leaflet(framenum,interval,nside,ntype,name_type,sside,nchain,carbons,ncarbons,array,odir):
  # sside : leafelt names
  # array : time series
  # odir  : output path
  for i in range(0,nside):
    side=sside[i]
    for j in range(0,ntype):
      ntail = nchain[j]
      for k in range(0,ntail):
        nc = ncarbons[j][k]

        # generate output string
        sout = f'# {side}: time series of {name_type[j]} SCD for chain {k}\n'
        sout +=f'#%10s' % carbons[j][k][0]
        for n in range(1,nc): sout+=f' {carbons[j][k][n]:10}'
        sout +=f'\n'

        for m in range(0,framenum):
          # sout += f'{interval*m}'
          for n in range(0,nc):
            sout += f' {array[i][j][k][n,m]:10.5f}'
          sout += f'\n'
        fout=f'{odir}/{side}_{name_type[j].lower()}_chain{k}.plo'
        f = open(fout,'w') ; f.write(sout); f.close() ;
        # print(sout)
  return



def write_time_series_bilayer(framenum,interval,ntype,name_type,nchain,carbons,ncarbons,array,odir):
  side="bilayer"
  for j in range(0,ntype):
    ntail = nchain[j]
    for k in range(0,ntail):
      nc = ncarbons[j][k]

      # generate output string
      sout = f'# {side}: time series of {name_type[j]} SCD for chain {k}\n'
      sout +='f#{carbons[j][k][0]:10s}'
      for n in range(1,nc): sout+=f' {carbons[j][k][n]:10s}' 
      sout +=f'\n'

      for m in range(0,framenum):
        # sout += f'{interval*m}'
        for n in range(0,nc):
          sout += f' {array[j][k][n,m]:10.5f}'
        sout += f'\n'
      fout=f'{odir}/{side}_{name_type[j].lower()}_chain{k}.plo'
      f = open(fout,'w') ; f.write(sout); f.close() ;
      # print(sout)
  return



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj','interval')
    parser.add_argument('--center', default=False, action='store_true',
      help='Perform membrane centering.')

    parser.add_argument('--sel', metavar='selection',
      help='Selection for individual molecule type with a format, "segid/resname/moltype MOLECULENAME and name CARBONNAMES or name REFATOMNAME". Reference atom is used for the leaflet assignment. For the CHARMM format topology, an selection for POPC can be "resname POPC and (name C22 or name C32 or name P)".')
    parser.add_argument('--split', action='store',# nargs='+',default='Y',
      help='Y/N. If N, the atom group for the selection is considered as that for a single molecule. If Y, the atom group is further splitted to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--sel-sys', metavar='selection',
      help='Selection for system atom groups in membranes for bilayer recentering.')
    parser.add_argument('--qz', action='store_true',default=False,
      help='Z-position based leaflet assigment: Default is false. Maybe useful when selections are minor components of the bilayer')
    parser.add_argument('--qb', action='store_true',default=False,
      help='Do bilayer analysis if provided: only for sym. bilayers')
    parser.add_argument('--qt', action='store_true',default=False,
      help='Write output time series if provided')
    parser.add_argument('--qa', action='store_true',default=False,
      help='Write output averages if provided')

    return parser



def run_scd(sel,split,qz,qb,qt,qa,sel_sys, psf:sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
  """
  ----------
  Calculate SCD order parameters.

  SCDs are calculated individual carbons in individual chains in individual molecules.
  Then, SCDs are avaraged for individual lipid types in individual leaflets.
  ----------
  """

  # process non-system arguments
  selection,ntype,sel_type,name_type,ref_name,ic_name,nchain,qsplit =\
    process_args(sel,split)

  # print summary
  for i in range(0,ntype):
    print(f'sel: {selection[i]} ; ref_name: {ref_name[i]}'\
          f' ; ic_name: {ic_name[i]} ; nchain = {nchain[i]}')
    print(f'#Split "{selection[i]}" into molecule level',qsplit[i])
  if (qz is True):
    print(f'Leaflets are assigned based on z-positions')
  else:
    print(f'Leaflets are assgined using a modified version of MDA LeafletFinder')
  print(f'Write average SCD,',qa) 
  print(f'Write time series of SCD',qt)
  print(f'SCD will be caulated every {interval} frames')
  
  # make output dir
  odir1="./scd/ave" ; os.system(f'mkdir -p {odir1}');
  odir2="./scd/time" ; os.system(f'mkdir -p {odir2}');
  
  # READ topology and trajectory
  u=mda.Universe(psf,traj) # MDA universe
  framenum=int(u.trajectory.n_frames/interval)  # number of frames to be analyzed

  # if (center is True): bilayer recentering - should be done before any assignments
  # - center in the box (an atom)
  # - center in the box (atom group for system)
  # - unwrap to get connectd molecules
  # Taken from 
  if  (center is True):
    origin = 0, 0, 0 #np.zeros([3],dtype=float) ; it did not work
    ag_cent = u.select_atoms(sel_sys)
    ag_all  = u.atoms
  
    workflow = [transformations.center_in_box(AtomGroup([ag_cent[0]]), point = origin),\
                transformations.center_in_box(ag_cent, point = origin),\
                transformations.unwrap(ag_all)]
  
    u.trajectory.add_transformations(*workflow)
  
  # Get numbers of molecules in individual lipid types,
  # Get total number of molecules,
  # Get lipid type index for individual molecules,
  # & Generate full atom groups for leaflet assignment
  #
  nmol_type,nmol,id_type,ag_full = \
    mymol.generate_full_mol_groups(u,ntype,sel_type,name_type,qsplit)
  
  # Generate
  #       carbon names for ind. chain in ind. lipid types
  #       hydrogen names boned to ind. carbons in ind. chain in ind. lipid types
  #       number of carbons in ind. chains in ind. lipid types
  #
  #	carbons[type][chain][carbon]
  #	hydrogens[type][chain][carbon][hydrogen]
  #	ncarbons[type][chain]
  #
  carbons,hydrogens,ncarbons = \
    myscd.generate_carbons_hydrogens_ncarbons_type(u,ag_full,ntype,name_type,\
                                                   nmol_type,ic_name,nchain)
  
  # Generate 
  # 	atomgroups of carbons in ind. chains in ind. lipids
  # 	ag_c[lipid][chain]
  #	ag_h[lipid][chain][carbon]
  ag_c,ag_h = \
    myscd.generate_ag_carbon_hydrogen(u,ag_full,name_type,id_type,\
                                      nchain,ncarbons,carbons,hydrogens)
  
  # Generate reference atom groups
  #	used to generate atom group for leaflet assignment
  #	used to assign leaflet index
  ag_ref = myscd.generate_ref_atomgroups(u,ag_full,ntype,nmol_type,ref_name)
  
  # raw SCD array for frames
  #  - raw data for individual carbon in indiviual tail in invidual lipids
  # raw_scd[lipid][chain][carbons] : list array upto lipid/chain level
  #                                : from carbons - numpy array
  #
  raw_scd = myscd.setup_raw_scd(nmol,id_type,nchain,ncarbons)
  
  # Set up output SCD
  scd,weight = myscd.setup_scd_weight(nside,ntype,nchain,ncarbons,framenum,qb) 
  ascd,astd = myscd.setup_ascd_astd(nside,ntype,nchain,ncarbons,qb) # average & std
  
  # Below, the membrane normal is assumed to be parallel to the z-axis
  memb_norm = np.array([0,0,1],dtype=float)
  
  # Analyze SCD
  # Loop over frames
  for i in range(0,framenum):
    print(f'# Analyze SCD: proceed {interval*i+1}/{interval*framenum}')
    ts = u.trajectory[interval*i]

    # do frame-wise bilayer recentering - remaining translation
    if (center is True):
      Lag_ref = myleaflet.assign_leaflet_zpos(u,ag_cent)
      zref = np.zeros([2],dtype=float)
      for i in range(0,nside):
        zref[i] = np.mean(Lag_ref[i].positions[:,2])
      # translation for z-centering
      tran = 0,0,-np.mean(zref)
      ts = transformations.translate(tran)(ts)
      ts = transformations.unwrap(ag_all)(ts)
  
    # Calculate Raw SCD
    myscd.calculate_raw_scd(raw_scd,nmol,ag_c,ag_h,memb_norm)
  
    # Assign leaflet indices for individual lipids: Run this even qb = False
    if (qz is True):
      if (i == 0): print(f'# leaflet assignemnt based on z-position')
      leaflets = myleaflet.assign_leaflet_zpos(u,ag_ref)
    else:
      if (i == 0): print(f'# leaflet assignemnt based on LeafletFinder with hydrid cutoff search')
      leaflets = myleaflet.assign_leaflet(u,ag_ref) 
    id_side = myscd.assign_leaflet_index(ag_ref,leaflets)
  
    # Now process to get scd for individual lipid types for the frame, i
    scd,weight = myscd.calculate_scd(i,ntype,raw_scd,scd,weight,id_side,id_type,nchain,ncarbons,qb)
  
  # Calculate weighted ave and std of SCD
  ascd,astd = myscd.calculate_ave_and_std(ascd,astd,scd,weight,framenum,qb)
  
  # Write output
  if (qb is True):
    if (qa is True):
      print(f'# Write average SCD')
      write_ave_std_bilayer(ntype,name_type,nchain,carbons,ncarbons,ascd,astd,odir1)
  
    if (qt is True):
      print(f'# Write time series output')
      write_time_series_bilayer(framenum,interval,ntype,name_type,nchain,carbons,ncarbons,scd,odir2)
  
  else: # if (qb is False):
    if (qa is True):
      print(f'# Write average SCD')
      write_ave_std_leaflet(nside,ntype,name_type,sside,nchain,carbons,ncarbons,ascd,astd,odir1)
  
    if (qt is True):
      print(f'# Write time series output')
      write_time_series_leaflet(framenum,interval,nside,ntype,name_type,sside,nchain,carbons,ncarbons,scd,odir2)
  
  return



def main(settings: Optional[dict] = None) -> None:
    if settings is None:
         settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_scd(**settings)



if __name__ == '__main__':
    main()

