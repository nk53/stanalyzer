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
from . import leaflet_util as myleaflet
from . import voronoi_analysis as myvorn
    
ANALYSIS_NAME = 'voronoi_shell_comp' 


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



def write_ave_std_leaflet(nside,ntype,nshell,sside,name_type,array1,array2,odir,otype):
  # sside: leaflet name
  # array1: average
  # array2: std
  # otype: output type name
  for i in range(0,nside):
    side = sside[i]
    for j in range(0,ntype):
      sout = f'# {side} {name_type[j]}: shell [{otype}_comp & std]_shell ntype= {ntype} nshell= {nshell}\n'
      sout+= f'#     shell'
      for m in range(0,ntype):
          sout += f' {name_type[m]:21s}'
      sout += f'\n'
      for k in range(1,nshell+1):
        sout += f' {k:10d}' # shell number
        for m in range(0,ntype):
          sout += f' {array1[i,j,k,m]:10.5f} {array2[i,j,k,m]:10.5f}'
        sout += f'\n'
      fout = f'{odir}/{side}_{name_type[j].lower()}_{otype}_comp.plo'
      f = open(fout,'w'); f.write(sout); f.close();
      print(sout)
  return



def write_ave_std_bilayer(ntype,nshell,name_type,array1,array2,odir,otype):
  # array1: average
  # array2: std
  # otype: output type name
  side = "bilayer"
  for j in range(0,ntype):
    sout = f'# {side} {name_type[j]}: shell [{otype}_comp & std]_shell ntype= {ntype} nshell= {nshell}\n'
    sout+= f'#     shell'
    for m in range(0,ntype):
        sout += f' {name_type[m]:21s}'
    sout += f'\n'
    for k in range(1,nshell+1):
      sout += f' {k:10d}' # shell number
      for m in range(0,ntype):
        sout += f' {array1[j,k,m]:10.5f} {array2[j,k,m]:10.5f}'
      sout += f'\n'
    fout = f'{odir}/{side}_{name_type[j].lower()}_{otype}_comp.plo'
    f = open(fout,'w'); f.write(sout); f.close();
    print(sout)
  return



def write_time_series_leaflet(framenum,interval,nside,ntype,nshell,sside,name_type,array,odir,otype):
  # sside: leaflet name
  # array: time series
  # otype: output type name
  for j in range(0,nside):
    side=sside[j]
    for k in range(0,ntype):
      sout = f'# {side} {name_type[k]}: time series of shell compositions\n'
      sout+= f'# time [{otype}_comp]_shell ntype= {ntype} nshell= {nshell}\n'
      sout+= f'#          '
      for m in range(1,nshell+1):
        tsshell = f'shell{m}'
        sout += f'{tsshell:44s}'
      sout += f'\n'
      sout += f'#          '
      for m in range(1,nshell+1):
        for n in range(0,ntype):
          sout += f' {name_type[n]:10s}'
      sout += f'\n'
      for i in range(0,framenum):
        # ct = (dt*(i+1)+(cnt-1))
        # sout += f' {ct}'
        sout += f' {interval*i:10d}'
        for m in range(1,nshell+1):
          for n in range(0,ntype):
            sout += f' {array[i,j,k,m,n]:10.5f}'
        sout += f'\n'
      fout = f'{odir}/time_{side}_{name_type[k].lower()}_{otype}_comp.plo'
      f=open(fout,'w'); f.write(sout); f.close();
      # print(sout)
  return



def write_time_series_bilayer(framenum,interval,ntype,nshell,name_type,array,odir,otype):
  # array: time series
  # otype: output type name
  side = "bilayer"
  for j in range(0,ntype):
    sout = f'# {side} {name_type[j]}: time series of shell compositions\n'
    sout+= f'# time [{otype}_comp]_shell ntype= {ntype} nshell= {nshell}\n'
    sout+= f'#          '
    for k in range(1,nshell+1):
      tsshell=f'shell{k}'
      sout += f'{tsshell:44s}'
    sout += f'\n'
    sout += f'#          '
    for k in range(1,nshell+1):
      for m in range(0,ntype):
        sout += f' {name_type[m]:10s}'
    sout += f'\n'
    for i in range(0,framenum):
      # ct = (dt*(i+1)+(cnt-1))
      # sout += f' {ct}'
      sout += f' {interval*i:10d}'
      for k in range(1,nshell+1):
        for m in range(0,ntype):
          sout += f' {array[i,j,k,m]:10.5f}'
      sout += f'\n'
    fout = f'{odir}/time_{side}_{name_type[j].lower()}_{otype}_comp.plo'
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



def run_voronoi_shell_comp(sel,split,qz,qb,qt,qa,sel_sys, psf:sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
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
  odir="./voronoi/shell";     os.system(f'mkdir -p {odir}');
  
  # READ topology and trajectory
  u=mda.Universe(psf,traj) # MDA universe
  framenum=int(u.trajectory.n_frames/interval) # number of frames to be analyzed
  
  # generate molecule atom groups
  name_type0,nmol_type0,nmol0,id_type0,ag0 = \
   mymol.generate_mol_groups(u,ntype,selection,qsplit)
  
  # atom group generation
  atomgroup = u.atoms[[]]
  for i in range(0,nmol0): atomgroup += ag0[i]
  
  ### Setup time series data
  # Shell info.: list array of time series
  # array[time][leaflet][type]...
  max_shell,t_numb_comp,t_frac_comp =\
    myvorn.setup_shell_comp_arrays(framenum,nside,ntype,qb)
  
  # START: LOOP over frames
  for i in range(0,framenum): 
    #  ct=(cnt-1)+dt*(i+1) # in ns
    print(f'# processing {interval*i+1}/{interval*framenum}')
    ts=u.trajectory[interval*i] 
  
    # get box size
    xtla,xtlb,xtlc=ts.dimensions[:3]
    box = np.array([xtla,xtlb],dtype=float)
  
    # initialize number of individual molecule types in the bilayer
    if (qb is True):
      nmol_typeb = np.zeros([ntype],dtype=float)
  
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
  
      # GENERATE Voronoi Center CONTACT INFO.
      print(f'# GENERATE VC CONTACT INFO.')
      contactVC,borderVC,distVC,weightVC =\
         myvorn.generate_vc_contact_info(vorn,nprim,nimage,box)
  
      # GENERATE MOLECULAR CONTACT INFO.
      print(f'# GENERATE MOL CONTACT INFO.')
      contactMol = \
        myvorn.generate_mol_contact_info(mol_ivc,contactVC,id_mol_vc,nprim)
      
      # SHELL-WISE NUMBER & FRACTIONAL COMPOSITIONS.
      print(f'# GENERATE SHELL-WISE NUMBER & FRACTIONAL COMPOSITIONS')
      for imol in range(0,nmol): 
        if (nmol >= 1000):
          if ( imol % int(nmol/10) == 0):
            print(f'shell analysis for molecule {imol}/{nmol}')
        itype = id_type[imol]    # molecule type index
        tnorm = nmol_type[itype] # number of the molecule type,itype
        tmax_shell,tnumb_comp,tfrac_comp = \
          myvorn.calculate_shell_comp(imol,nmol,contactMol,ntype,mol_ivc,id_type) 
  
        # update time series (un-normalized)
        # Accumulating numbers.
        if (qb is True):
          max_shell[i][itype],t_numb_comp[i][itype],t_frac_comp[i][itype] = \
            myvorn.update_shell_comp(ntype,max_shell[i][itype],\
                          t_numb_comp[i][itype],t_frac_comp[i][itype],\
                          tmax_shell,tnumb_comp,tfrac_comp)
        else: # if (qb is False):
          max_shell[i][iside][itype],\
          t_numb_comp[i][iside][itype],t_frac_comp[i][iside][itype] = \
            myvorn.update_shell_comp(ntype,max_shell[i][iside][itype],\
                          t_numb_comp[i][iside][itype],t_frac_comp[i][iside][itype],\
                          tmax_shell,tnumb_comp,tfrac_comp)
  
      if (qb is True):
        # update nmol_typeb for normalization
        nmol_typeb += nmol_type
      else:
      # if (qb is False):
        t_numb_comp[i][iside],t_frac_comp[i][iside] = \
          myvorn.normalize_raw_comp(t_numb_comp[i][iside],t_frac_comp[i][iside],nmol_type)
    # normalize bilayer data
    if (qb is True):
      t_numb_comp[i],t_frac_comp[i] = \
        myvorn.normalize_raw_comp(t_numb_comp[i],t_frac_comp[i],nmol_typeb)
  # END: LOOP over frames
  
  # STAT for average shell compositions
  # set a shell number, nshell
  # results will be saved upto the shell, nshell
  nshell = myvorn.set_nshell(max_shell,framenum,nside,ntype,qb)
  
  # Process to get statistics...
  if (qa is True):
    # Convert time series list arrays to ndarrays
    numb_comp,frac_comp = myvorn.convert_to_ndarray(t_numb_comp,t_frac_comp,framenum,nside,ntype,nshell,qb)
  
    anumb_comp = np.average(numb_comp,axis=0)
    snumb_comp = np.std(numb_comp,axis=0)
    
    afrac_comp = np.average(frac_comp,axis=0)
    sfrac_comp = np.std(frac_comp,axis=0)
  
    print(f'# Write average output')   
    if (qb is True):
      # write NUMB_COMP outputs # 1 to nshell
      otype="numb"
      write_ave_std_bilayer(ntype,nshell,name_type,anumb_comp,snumb_comp,odir,otype)
  
      # write FRAC_COMP outputs: shell 1 to nshell
      otype="frac"
      write_ave_std_bilayer(ntype,nshell,name_type,afrac_comp,sfrac_comp,odir,otype)
    else:
    # if (qb is False):
      # write NUMB_COMP outputs # 1 to nshell
      otype="numb"
      write_ave_std_leaflet(nside,ntype,nshell,sside,name_type,anumb_comp,snumb_comp,odir,otype)
      
      # write FRAC_COMP outputs: shell 1 to nshell
      otype="frac"
      write_ave_std_leaflet(nside,ntype,nshell,sside,name_type,afrac_comp,sfrac_comp,odir,otype)
  
  # Time series output 
  if (qt is True): # pass
    print(f'# Write time series output')
    # dt=1.0/float(framenum) # increment in time
    if (qb is True):
      otype="numb"
      write_time_series_bilayer(framenum,interval,ntype,nshell,name_type,numb_comp,odir,otype)
      
      # fractional composition/shell
      otype="frac"
      write_time_series_bilayer(framenum,interval,ntype,nshell,name_type,frac_comp,odir,otype)
    else:
    # if (qb is False):
      otype="numb"
      write_time_series_leaflet(framenum,interval,nside,ntype,nshell,sside,name_type,numb_comp,odir,otype)
    
      # fractional composition/shell
      otype="frac"
      write_time_series_leaflet(framenum,interval,nside,ntype,nshell,sside,name_type,frac_comp,odir,otype)

  return



def main(settings: Optional[dict] = None) -> None:
    if settings is None:
         settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_voronoi_shell_comp(**settings)



if __name__ == '__main__':
    main()


