#!/usr/bin/python
import argparse
from typing import Optional
import os,sys,re,math
import stanalyzer.bin.stanalyzer as sta
import numpy as np                              
import MDAnalysis as mda  
import MDAnalysis.transformations as transformations
from MDAnalysis.core.groups import AtomGroup
from . import mol_atomgroup_util as mymol
from . import msd_util as mymsd

ANALYSIS_NAME = 'msd_membrane'


#--- The following are hard set for membrane analysis
nside=2     # up/dn
sside=["up","dn"]



def process_args(sel,split,sel_sys):
  """
  ----------
  Process arguments
  ----------
  """

  selection = re.split(';|,', f'{sel:s}')
  ntype = len(selection)
  for i in range(0,ntype): selection[i]=selection[i].strip()
  split = re.split(';|,',f'{split:s}')
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

  sel_sys = f'{sel_sys:s}'
  return(selection,ntype,qsplit,sel_sys)



# Write system COM
def write_com_sys(traj_com_sys_unwrap,nside,framenum,interval,odir,sside):
  for i in range(0,nside):
    sout =f'#     frame       COM\n'
    sout+=f'#                   x           y          z\n'
    for j in range(0,framenum):
      tcom = traj_com_sys_unwrap[i][j]
      sout+=f' {interval*j:10d} {tcom[0]:10.5f} {tcom[1]:10.5f} {tcom[2]:10.5f}\n'
    # print(sout)
    fout=f'{odir}/{sside[i]}_com_sys_unwrapped.dat'
    f=open(fout,'w'); f.write(sout); f.close();
  return



# Write unwrapped COMs of individual molecule
def write_com_mol(traj_com_unwrap,nside,nmol,framenum,interval,odir,sside):
  for i in range(0,nside):
    sout =f'#  leaflet {sside[i]}\n'
    sout+=f'#     frame       COM\n'
    sout+=f'#                   x           y          z\n'
    for j in range(0,framenum):
      tcom = traj_com_unwrap[i][j]
      sout+=f' {interval*j:10d}' 
      for k in range(0,nmol[i]):
        for m in range(0,3):
          sout+=f' {tcom[k,m]:10.5f}' 
      sout+=f'\n'
    # print(sout)
    fout=f'{odir}/{sside[i]}_com_mol_unwrapped.dat'
    f=open(fout,'w'); f.write(sout); f.close();
  return



# Write x,y,z-components of MSD for given molecule type in a given leaflet
def write_msd(name,msd,taus,odir,side):
  sout = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'

  ntau = len(taus)
  for i in range(0,ntau):
    sout += f' {taus[i]:10.5f}'
    for j in range(0,3):
      sout += f' {msd[i,j]:10.5f}'
    sout += '\n'

  fout = f'{odir}/{side}_{name.lower()}_msd.dat'
  f = open(fout,'w') ; f.write(sout); f.close();
  return



# Write MSD outputs for leaflets
def write_msd_outputs_leaflet(msd,taus,nside,ntype,name_type,odir,sside):
  for i in range(0,nside):
    side = sside[i]

    for j in range(0,ntype):
      name = name_type[j]

      print(f'# Write MSDs for {name} in leafelt, {side}')
      write_msd(name,msd[i][j],taus,odir,side)

  return



# Write MSD outputs for bilayer
def write_msd_outputs_bilayer(bmsd,taus,ntype,name_type,odir):
  side = "bilayer"

  for i in range(0,ntype):
    name = name_type[i]

    print(f'# Write MSDs for {name} in bilayer')
    write_msd(name,bmsd[i],taus,odir,side)

  return



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj','interval')
    parser.add_argument('--center', default=False, action='store_true',
      help='Perform membrane centering.')

    parser.add_argument('--sel', metavar='selection', # nargs='+',
      help='Selection for individual molecule. MSD will be calculated for COMs of individual molecules.')
    parser.add_argument('--split', action='store',# nargs='+',default='Y',
      help='Y/N. If N, the atom group for the selection is considered as that for a single molecule. If Y, the atom group is further splitted to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--sel-sys', metavar='selection',
      help='Selection for system atom groups for leaflet COM drift correction and bilayer recentering.')
    parser.add_argument('--qz', action='store_true',default=False,
      help='Z-position based leaflet assigment: Default is false. Maybe useful when selections are minor components of the bilayer')
    parser.add_argument('--qb', action='store_true',default=False,
      help='Do bilayer analysis if provided: only for sym. bilayers')
    parser.add_argument('--qcomsys',default=False,action='store_true',
      help='When provided, write unwrapped COM trajectories for leaflets.')
    parser.add_argument('--qcommol',default=False,action='store_true',
      help='When provided, write unwrapped COM trajectoreis for molecules.')

    return parser



def run_msd_membrane(sel,split,sel_sys,qz,qb,qcomsys,qcommol, psf:sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
  """
  ----------
  Calculate Mean Square Displacement of COMs of selected molecule types

  MSDs are calculated separately for individual leaflets.
  Results will be obtained for leaflets or bilayer (weighted average of leaflet MSDs).
  ----------
  """

  # process arguments
  selection,ntype,qsplit,sel_sys = process_args(sel,split,sel_sys)


  # print summary of arguments
  for i in range(0,ntype):
    print(f'#Split "{selection[i]}" into molecule level',qsplit[i])
  if (qz is True):
    print(f'Leaflets are assigned based on z-positions')
  else:
    print(f'Leaflets are assgined using a modified version of MDA LeafletFinder')
  print(f'Writing unwrapped COM of leaflets:',qcomsys)
  print(f'Writing unwrapped COM of individual molecules:',qcommol)
  print(f'MSD will be caulated every {interval} frames in delay time')
  print(f'Bilayer is recentered at z = 0 using {sel_sys}:',center)
  
  odir = "msd" ; os.system(f'mkdir -p {odir}');
  
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
  
  # LEAFLETs are assigned in this stage and will not be altered.
  # Generate molecule groups
  name_type,nmol_type,nmol,id_type,ag =\
    mymol.generate_mol_groups_memb(u,nside,ntype,selection,qsplit,sside,qz)
  
  # Generate system groups: required for leaflet COM drift correction
  ag_sys = mymol.generate_sys_groups(u,sel_sys,qz)
  
  # Set arrays in use
  pos,pos_prev = [],[]                # current and previous atom positions in ind. molecules
  displ,pos_unwrap = [],[]            # displacements and unwrapped position of atoms
  mass_mol,tmass_mol = [],[]          # atom masses and total mass of ind. molecules
  com_unwrap,traj_com_unwrap = [],[]  # unwrapped COM and its traj of ind. molecules
  
  pos_sys,pos_sys_prev = [],[]        # current and previous atom positions in ind. leaflets
  mass_sys,tmass_sys = [],[]          # atom masses and total mass of ind. leaflets
  displ_sys = []                      # displacements of atoms in individual leaflets
  displ_sys_com = []                  # COM displacement of individual leaflets
  com_sys_unwrap = []                 # unwrapped leaflet COM
  traj_com_sys_unwrap = []            # trajectory of unwrapped leaflet COM
  
  for i in range(0,nside):
    # Set mass, position, and displacement arrays 
    # Will be input for calculation of unwrapped COM of individual molecules
    mass_mol.append([]) ; tmass_mol.append([]);
    pos.append([]);       pos_prev.append([]);
    displ.append([]);     pos_unwrap.append([]);
    mass_mol[i],tmass_mol[i],pos[i],pos_prev[i],displ[i],pos_unwrap[i] =\
      mymsd.set_mass_pos_displ_arrays(nmol[i],ag[i])
  
    # Set unwrapped COM positions and their trajectories for individual molecules
    com_unwrap.append([]); traj_com_unwrap.append([]);
    com_unwrap[i],traj_com_unwrap[i] =\
      mymsd.setup_unwrapped_mol_com_traj_array(ag[i],framenum)
  
    # set system masses, positions, dsplacements
    mass_sys.append([]);  tmass_sys.append([]);
    pos_sys.append([]);   pos_sys_prev.append([]);
    displ_sys.append([]); displ_sys_com.append([]);
  
    mass_sys[i],tmass_sys[i],pos_sys[i],pos_sys_prev[i],\
    displ_sys[i],displ_sys_com[i] = \
      mymsd.setup_sys_mass_pos_displ_arrays(ag_sys[i])
    
    # set unwrapped system COM arrays
    com_sys_unwrap.append([]); traj_com_sys_unwrap.append([]);
    com_sys_unwrap[i],traj_com_sys_unwrap[i] =\
      mymsd.setup_unwrapped_com_traj_array(framenum)
  
  # UNWRAPPING
  print(f'# UNWRAP trajectories')
  ## sys.exit(0)
  for i in range(0,framenum):
    #  ct=(cnt-1)+dt*(i+1) # in ns
    print(f'# processing {interval*i+1}/{interval*framenum}')
    ts=u.trajectory[interval*i] # Is this update frame ?

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
  
    # get box size
    xtla,xtlb,xtlc=ts.dimensions[:3]
    box = np.array([xtla,xtlb,xtlc],dtype=float)
  
    # read cooridnates 
    for j in range(0,nside):
      mymsd.read_coor(nmol[j],pos[j],ag[j])
      pos_sys[j] = ag_sys[j].positions
  
      if (i==0):
        ## Initialization ##
        # positions, unwrapped COM, and traj of unwrapped COM for the entire system 
        mymsd.init_unwrap_sys_com(pos_sys[j],mass_sys[j],tmass_sys[j],\
                            pos_sys_prev[j],com_sys_unwrap[j],traj_com_sys_unwrap[j])
        
        # positions, unwrapped COM, and traj of unwraped COM for individual molecules
        mymsd.init_unwrap_mol_com(pos[j],mass_mol[j],tmass_mol[j],pos_prev[j],\
                            pos_unwrap[j],com_unwrap[j],traj_com_unwrap[j])
        # print(f'# leaflet {sside[j]}: init unwrapped com/traj done')
      else:
        # get sys COM displacement and update positions of system (in pos_sys_prev)
        # Current positions of system is obtained inside the function
        mymsd.calculate_displ_sys_com(i,box,pos_sys[j],pos_sys_prev[j],\
                                      mass_sys[j],tmass_sys[j],displ_sys_com[j])
        # print(f'# leafelt {sside[j]}: sys COM displacement:',displ_sys_com[j])
  
        # update unwrapped system COM 
        com_sys_unwrap[j] = com_sys_unwrap[j] + displ_sys_com[j]
        # update unwrapped system COM trajectory
        np.copyto(traj_com_sys_unwrap[j][i],com_sys_unwrap[j])
  
        # update unwrapped mol positions
        mymsd.update_unwrapped_mol_pos(i,box,pos[j],pos_prev[j],\
                                       pos_unwrap[j],displ_sys_com[j])
  
        # calculate mol COMs & update their traj
        mymsd.update_unwrapped_mol_com_traj(i,pos_unwrap[j],mass_mol[j],tmass_mol[j],\
                                            com_unwrap[j],traj_com_unwrap[j])
  
        # print(com_sys_unwrap[j])
        # print(com_unwrap[j])
  
  print(f'# UNWRAPPING TRAJ & COMS DONE')
  # sys.exit(0)
  
  if (qcomsys is True):
    print(f'# Write unwrapped COM of leaflets')
    write_com_sys(traj_com_sys_unwrap,nside,framenum,interval,odir,sside)
  
  if (qcommol is True):
    print(f'# Write unwrapped COMs of individual molecules')
    write_com_mol(traj_com_unwrap,nside,nmol,framenum,interval,odir,sside)
  
  print(f'# MSD calculations')
  # Loop over delay times with given interval
  taus = []; [taus.append(interval*i) for i in range(0,framenum)];
  ntau = len(taus) # number of data points along the delay time
  
  # Setup msd for individual molecule types
  msd = []
  for i in range(0,nside):
    tmsd = mymsd.setup_msd_arrays(ntype,ntau)
    msd.append(tmsd)
  
  # Calculate MSD for delay times, tau in {taus}
  for i in range(0,nside):
    print(f'# leaflet {sside[i]}')
    mymsd.calculate_msd(taus,framenum,traj_com_unwrap[i],id_type[i],msd[i])
  
  # Write MSD outputs
  if (qb is True):
    # calculate MSD over bilayers
    bmsd = mymsd.calculate_msd_bilayer(msd,nside,ntype,ntau,nmol_type)
    write_msd_outputs_bilayer(bmsd,taus,ntype,name_type,odir)
  else: # if (qb is False):
    write_msd_outputs_leaflet(msd,taus,nside,ntype,name_type,odir,sside)
  
  return



def main(settings: Optional[dict] = None) -> None:
    if settings is None:
         settings = dict(sta.get_settings(ANALYSIS_NAME))
    # non-system arguments will be handled at the beginnig of this function
    run_msd_membrane(**settings)



if __name__ == '__main__':
    main()


