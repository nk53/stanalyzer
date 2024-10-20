#!/usr/bin/python
import argparse
from typing import Optional
import os,sys,re,math
import stanalyzer.bin.stanalyzer as sta
import numpy as np                                
import MDAnalysis as mda                          
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis.analysis.leaflet as mdaleaflet  
from scipy.spatial import Voronoi                   
from . import mol_atomgroup_util as mymol
from . import voronoi_analysis as myvorn

ANALYSIS_NAME = 'voronoi_apl' 


#--- The following are hard set for membrane analysis
dim=2              # 2 dimension
nimage=int(3**dim) # number of total primary+images for 2D

nside=2     # up/dn
sside=["up","dn"]



def process_args(sel,split):
  """
  ---------- 
  Process arguments
  ----------
  """ 

  selection = re.split(';|,', f'{sel:s}') # selections for individual molecule types 
  ntype = len(selection)                  # number of unique molecule types
  for i in range(0,ntype): selection[i]=selection[i].strip()

  split  = re.split(';|,',f'{split:s}')
  nsplit = len(split)
  for i in range(0,nsplit): split[i]=split[i].strip()
  qsplit = []
  if (nsplit < ntype): # add more qsplit options
    for i in range(0,nsplit):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)
    for i in range(nsplit,ntype): qsplit.append(True) # default values
  else: #  get split upto ntype
    qsplit = []
    for i in range(0,ntype):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)

  return(selection,ntype,qsplit)



def write_ave_std_leaflet(nside,ntype,name_type,sside,array1,array2,array3,odir):
  # array1: average
  # array2: std
  # array3: total weight
  for i in range(0,nside):
    side=sside[i]
    sout = f'# {side}: APL & std\n#'
    for j in range(0,ntype):
      sout += f' {name_type[j]:21s}'
    sout += f'\n'
    for j in range(0,ntype):
      sout += f' {array1[i,j]:10.5f} {array2[i,j]:10.5f}'
    sout += f'\n'
    # add numbers
    sout += f'# weights:'
    for j in range(0,ntype):
      sout += f' {array3[i,j]:10.5f}'

    fout = f'{odir}/{side}.plo'
    f=open(fout,'w') ; f.write(sout); f.close();
    print(sout)
  return



def write_ave_std_bilayer(ntype,name_type,array1,array2,array3,odir):
  # array1: average
  # array2: std
  # array3: total weight
  side="bilayer"
  sout = f'# {side}: APL & std\n#'
  for j in range(0,ntype):
    sout += f' {name_type[j]:21s}'
  sout += f'\n'
  for j in range(0,ntype):
    sout += f' {array1[j]:10.5f} {array2[j]:10.5f}'
  sout += f'\n'
  # add numbers
  sout += f'# weights:'
  for j in range(0,ntype):
    sout += f' {array3[j]:10.5f}'

  fout = f'{odir}/{side}.plo'
  f=open(fout,'w') ; f.write(sout); f.close();
  print(sout)
  return



def write_time_series_leaflet(framenum,interval,nside,ntype,name_type,sside,array,odir):
  # sside: leaflet names
  # array: time series
  # odir : output path
  for j in range(0,nside):
    side=sside[j]
    sout  = f'# {side}: APL time series\n'
    sout += f'#     frame'
    for k in range(0,ntype):
      sout += f' {name_type[k]:10s}'
    sout += f'\n'
    for i in range(0,framenum):
      # ct = (dt*(i+1) + (cnt - 1)) # 
      # sout+= f' {ct}'
      sout += f' {interval*i:10d}'
      for k in range(0,ntype):
        sout += f' {array[i,j,k]:10.5f}'
      sout += f'\n'

    fout = f'{odir}/time_{side}.plo'
    f=open(fout,'w'); f.write(sout); f.close();
    # print(sout)
  return



def write_time_series_bilayer(framenum,interval,ntype,name_type,array,odir):
  # array: time series
  # odir : output path
  side="bilayer"
  sout = f'# {side}: APL time series\n'
  sout+= f'#     frame'
  for j in range(0,ntype):
    sout += f' {name_type[j]:10s}'
  sout += f'\n'
  for i in range(0,framenum):
    # ct = (dt*(i+1) + (cnt - 1)) # 
    # sout += f' {ct}'
    sout += f' {interval*i:10d}'
    for j in range(0,ntype):
      sout += f' {array[i,j]:10.5f}'
    sout += f'\n'

  fout = f'{odir}/time_{side}.plo'
  f=open(fout,'w'); f.write(sout); f.close();
  # print(sout)
  return


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj','interval')
    parser.add_argument('--center', default=False, action='store_true',
      help='Perform membrane centering.')

    parser.add_argument('--sel', metavar='selection',
      help='Selection for individual molecule type with a format, "segid/resname/moltype MOLECULENAME and name ATOMNAMES". For the CHARMM format topology, an selection for POPC can be "resname POPC and (name C2 or name C21 or name C31)".')
    parser.add_argument('--split', action='store',
      help='Y/N. If N, the atom group for the selection is considered as that for a single molecule. If Y, the atom group is further splitted to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--sel_sys', metavar='selection',
      help='Selection for system atom groups for membrane recentering.')
    parser.add_argument('--qz', action='store_true',default=False,
      help='Z-position based leaflet assigment: Default is false. Maybe useful when selections are minor components of the bilayer')
    parser.add_argument('--qb', action='store_true',default=False,
      help='Do bilayer analysis if provided: only for sym. bilayers')
    parser.add_argument('--qt', action='store_true',default=False,
      help='Write output time series if provided')
    parser.add_argument('--qa', action='store_true',default=False,
      help='Write output averages if provided')

    return parser



def run_voronoi_apl(sel,split,qz,qb,qt,qa,sel_sys, psf:sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
  """
  ----------
  Calculate area per lipid (APL) of inidividual molecul types.

  APLs are calculated individual molecules in individual leaflets.
  Then, APLs are avaraged for individual lipid types.
  ----------
  """

  # process non-system arguments
  selection,ntype,qsplit = process_args(sel,split)

  # print summary of arguments
  for i in range(0,ntype):
    print(f'#Split "{selection[i]}" into molecule level',qsplit[i])
  print(f'Write results for bilayer',qb)
  print(f'Write averge APLs',qa)
  print(f'Write time sereis output',qt)
  print(f'Bilayer is recentered at z = 0 using {sel_sys}:',center)
  print(f'Analysis will be done every {interval}frames')
  
  # make output dir
  odir="./voronoi/apl";       os.system(f'mkdir -p {odir}');
  
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
  
  # Generate molecule atom groups
  name_type0,nmol_type0,nmol0,id_type0,ag0 = \
   mymol.generate_mol_groups(u,ntype,selection,qsplit)
  
  # Atom group for all atoms
  atomgroup = u.atoms[[]]
  for i in range(0,nmol0):
    atomgroup += ag0[i] 
  
  ### Setup time series data
  if (qb is True):
    apl = np.zeros([framenum,ntype],dtype=float)
    ncomp = np.zeros([framenum,ntype],dtype=float)
  else: # if (qb is False):
    apl = np.zeros([framenum,nside,ntype],dtype=float)
    ncomp = np.zeros([framenum,nside,ntype],dtype=float)
  
  # START: LOOP over frames
  for i in range(0,framenum): 
    #  ct=(cnt-1)+dt*(i+1) # in ns
    print(f'# processing {interval*i+1}/{interval*framenum}')
    ts=u.trajectory[interval*i] 

    # do frame-wise bilayer recentering - remaining translation
    if (center is True):
      Lag_ref = myleaflet.assign_leaflet(u,ag_cent)
      zref = np.zeros([2],dtype=float)
      for i in range(0,nside):
        zref[i] = np.mean(Lag_ref[i].positions[:,2])
      # translation for z-centering
      tran = 0,0,-np.mean(zref)
      ts = transformations.translate(tran)(ts)
      ts = transformations.unwrap(ag_all)(ts)
  
    # get box size
    xtla,xtlb,xtlc=ts.dimensions[:3]
    box = np.array([xtla,xtlb],dtype=float)
  
    # Assign leaflet & generate mol groups for membrane
    name_type,Lnmol_type,Lnmol,Lid_type,Lag = \
      mymol.generate_mol_groups_memb(u,nside,ntype,selection,qsplit,sside,qz)
  
    # START: LOOP over leaflets
    for iside in range(0,nside):
      print(f'# leaflet {iside}')
  
      nmol = Lnmol[iside]             # number of molecules in the leaflet
      nmol_type = Lnmol_type[iside]   # nmuber of molecules for individual types 
      id_type = Lid_type[iside]       # molecule type index for individual molecules
  
      # generate VC list for individual molecules &
      # assign molecule index to individual VCs
      mol_ivc,id_mol_vc =\
         myvorn.connect_mol_and_vc(u,nmol,nmol_type,Lag[iside])
  
      # PRINT BASIC INFO.
      sinfo  = f'# SYS. INFO.\n'
      sinfo += f'\t numb. of user defined mol. types= {ntype}\t'
      for itype in range(0,ntype): sinfo += f' {name_type[itype]}'
      sinfo += f'\n\t\t\t\t\t\t'
      for itype in range(0,ntype): sinfo += f' {nmol_type[itype]:4d}'
      sinfo += f' ntot= {nmol}\n'
      sinfo += f'\t system size\t\t\t\t xtla= {xtla:10.5f} xtlb= {xtlb:10.5f}'
      print(sinfo)
  
      # PREPARE VORONOI CENTER
      print(f'# VORONOI TESSELLATION')
      nprim,nvorn,vc = \
        myvorn.prepare_voronoi_centers(Lag[iside],box,dim,nimage)
    
      # VORONOI TESSEELATION
      vorn=Voronoi(vc)
  
      # CALCULATE AREA OF EACH MOLECULE
      area = myvorn.calc_mol_area_all(nmol,mol_ivc,vorn)
  
      # CALCULATE COMPONENT APL
      tapl,tncomp = myvorn.calc_apl(area,id_type,ntype) 
      # update time series data
      if (qb is True):
        apl[i,:] += tncomp * tapl
        ncomp[i,:] += tncomp
      else:  # if (qb is False):
        np.copyto(apl[i,iside,:],tapl)
        np.copyto(ncomp[i,iside,:],tncomp)
  
  # END: LOOP over frames
  if (qb is True):
    # get APL 
    for i in range(0,framenum):
      apl[i,:] /= ncomp[i,:]
  
  # Process to get statistics...
  if (qa is True):
    # aapl: average component APL
    # sapl: std of component APL
    # anum: total counts for component APL
    aapl,sapl,anum = myvorn.calculate_ave_and_std(nside,ntype,apl,ncomp,qb)
    
    # write APL outputs
    print(f'# Write average ouput')
    if (qb is True):
      write_ave_std_bilayer(ntype,name_type,aapl,sapl,anum,odir)
    else: # if (qb is False):
      write_ave_std_leaflet(nside,ntype,name_type,sside,aapl,sapl,anum,odir)
       
  # Time series output 
  if (qt is True): # pass
    print(f'# Write time series output')
    # dt=1.0/float(framenum) # increment in time
    if (qb is True):
      write_time_series_bilayer(framenum,interval,ntype,name_type,apl,odir)
    else: # if (qb is False):
      write_time_series_leaflet(framenum,interval,nside,ntype,name_type,sside,apl,odir)
  
  return



def main(settings: Optional[dict] = None) -> None:
    if settings is None:
         settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_voronoi_apl(**settings)



if __name__ == '__main__':
    main()


