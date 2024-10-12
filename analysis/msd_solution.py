#!/usr/bin/python
import argparse
from typing import Optional
import os,sys,re
import stanalyzer.bin.stanalyzer as sta
import numpy as np                 
import MDAnalysis as mda           
from . import mol_atomgroup_util as mymol
from . import msd_util as mymsd

ANALYSIS_NAME = 'msd_solution'



def process_args(sel,split):
  """
  ----------
  Process arguements
  ----------
  """

  selection = re.split(';|,', f'{sel:s}')
  ntype = len(selection)
  for i in range(0,ntype): selection[i]=selection[i].strip()
  # print("selection",selection);print(ntype); sys.exit(0);
  split = re.split(';|,',f'{split:s}')
  nsplit = len(split)
  for i in range(0,nsplit): split[i]=split[i].strip()
  qsplit = []
  if (nsplit < ntype): # add more qsplit options
    for i in range(0,nsplit):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)
    for i in range(nsplit,ntype): qsplit.append(True) # default values
  else: #  get my_dict["split"] upto ntype
    qsplit = []
    for i in range(0,ntype):
      if (split[i].lower() == "y"): qsplit.append(True)
      else: qsplit.append(False)
  return(selection,ntype,qsplit)  



def write_com_mol(traj_com_unwrap,nmol,framenum,odir):
  """
  ----------
  Write unwrapped COMs of individual molecules
  ----------
  """

  sout =f'#     frame       COM\n'
  sout+=f'#                   x           y          z\n'
  for i in range(0,framenum):
    tcom = traj_com_unwrap[i]
    sout+=f' {i:10d}'
    for j in range(0,nmol):
      for k in range(0,3):
        sout+=f' {tcom[j,k]:10.5f}'
    sout+=f'\n'
  # print(sout)
  fout=f'{odir}/com_mol_unwrapped.dat'
  f=open(fout,'w'); f.write(sout); f.close();
  return



# Write x,y, and z-components of MSD for given molecule type
def write_msd(name,msd,taus,odir):
  """
  ----------
  Write x,y, & z-components of MSD for a given molecule type"
  ----------
  """

  sout = f'#{"tau":10s} {"MSDX":10s} {"MSDY":10s} {"MSDZ":10s}\n'
  ntau = len(taus)
  for i in range(0,ntau):
    sout += f'{taus[i]:10.5f}'
    for j in range(0,3): # x,y,z
      sout += f' {msd[i,j]:10.5f}'
    sout += f'\n'
  
  fout = f'{odir}/{name.lower()}_msd.dat'
  f = open(fout,'w'); f.write(sout); f.close();
  return 

# Write MSD outputs
def write_msd_outputs(msd,taus,ntype,name_type,odir):
  for i in range(0,ntype):
    name = name_type[i]
    print(f'# Write MSD output for {name}')
    write_msd(name,msd[i],taus,odir)

  return



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj','interval')
    parser.add_argument('--center', default=False, action='store_true',
      help='Perform membrane centering.')

    parser.add_argument('--sel', metavar='selection',# nargs='+',
      help='Selection for individual molecule. MSD will be calculated for COMs of individual molecules.')
    parser.add_argument('--split', action='store',# nargs='+',default='Y',
      help='Y/N. If N, the atom group for the selection is considered as that for a single molecule. If Y, the atom group is further splitted to molecular level based on segid/resname/moleculename. Default is Y/y.')
    parser.add_argument('--qcom', default=False, action='store_true',
      help='When provided, write a COMs trajectory for molecules.')

    return parser


def run_msd_solution(sel,split,qcom, psf:sta.FileRef, traj: sta.FileRefList, interval: int = 1, center: bool = False) -> None:
  """
  ----------
  Calculate Mean Square Displacement of COMs of selected molecule types
  ----------
  """

  # process non-system arguments
  selection,ntype,qsplit = process_args(sel,split)

  # print summary of arguments
  for i in range(0,ntype):
    print(f'#Split "{selection[i]}" into molecule level',qsplit[i])
  print(f'Writing unwrapped COM of individual molecules:',qcom)
  print(f'MSD will be caulated every {interval} frames in delay time')
  
  odir = "msd" ; os.system(f'mkdir -p {odir}');
  
  # READ topology and trajectory
  u=mda.Universe(psf,traj) # MDA universe
  framenum=u.trajectory.n_frames  # number of frames
  
  # Generate molecule groups
  name_type,nmol_type,nmol,id_type,ag =\
    mymol.generate_mol_groups(u,ntype,selection,qsplit)
  
  # Set mass, position, and displacement arrays 
  # Will be input for calculation of unwrapped COM of individual molecules
  mass_mol,tmass_mol,pos,pos_prev,displ,pos_unwrap =\
    mymsd.set_mass_pos_displ_arrays(nmol,ag)
  
  # Set unwrapped COM positions and their trajectories for individual molecules
  com_unwrap,traj_com_unwrap =\
    mymsd.setup_unwrapped_mol_com_traj_array(ag,framenum)
  
  ## For solution system (There should be no system COM drift)
  displ_sys_com = np.zeros([3],dtype=float)
  
  ## ag_sys = u.atoms[[]]
  ## for i in range(0,nmol): ag_sys += ag[i]
  ## 
  ## mass_sys,tmass_sys,pos_sys,pos_sys_prev,displ_sys,displ_sys_com = \
  ##   mymsd.setup_sys_mass_pos_displ_arrays(ag_sys)
  ## 
  ## com_sys_unwrap,traj_com_sys_unwrap = \
  ##   mymsd.setup_unwrapped_com_traj_array(ag_sys,framenum)
  
  # UNWRAPPING
  print(f'# UNWRAP trajectories')
  for i in range(0,framenum):
    #  ct=(cnt-1)+dt*(i+1) # in ns
    print(f'# processing {i+1}/{framenum}')
    ts=u.trajectory[i] # Is this update frame ?
  
    # get box size
    xtla,xtlb,xtlc=ts.dimensions[:3]
    box = np.array([xtla,xtlb,xtlc],dtype=float)
  
    # read cooridnates (translate by ?)
    mymsd.read_coor(nmol,pos,ag)
  
    if (i==0):
      ## Initialization ##
      ## mymsd.init_unwrap_sys_com(pos_sys,mass_sys,tmass_sys,pos_sys_prev,\
      ##                           com_sys_unwrap,traj_com_sys_unwrap)
  
      # positions, unwrapped COM, and traj of unwraped COM for individual molecules
      mymsd.init_unwrap_mol_com(pos,mass_mol,tmass_mol,pos_prev,\
                                pos_unwrap,com_unwrap,traj_com_unwrap)
      print(f'# init unwrap_pos done')
      # sys.exit(0)
    else:
      ## # TESTED - SYS COM show no displacement
      ## # calculate sys COM displacement
      ## mymsd.calculate_displ_sys_com(i,box,pos_sys,pos_sys_prev,mass_sys,tmass_sys,\
      ##                               displ_sys_com)
      ## print(f'sys COM displacement:',displ_sys_com)
      ## # update unwrapped system COM
      ## tcom_sys_unwrap = com_sys_unwrap + displ_sys_com
      ## np.copyto(com_sys_unwrap,tcom_sys_unwrap)
      ## # update unwrapped system COM trajectory
      ## np.copyto(traj_com_sys_unwrap[i],com_sys_unwrap)
  
      # update unwrapped mol positions
      mymsd.update_unwrapped_mol_pos(i,box,pos,pos_prev,pos_unwrap,displ_sys_com)
  
      # calculate mol COMs & update their traj
      mymsd.update_unwrapped_mol_com_traj(i,pos_unwrap,mass_mol,tmass_mol,com_unwrap,traj_com_unwrap)
  
  print(f'# UNWRAPPING TRAJ & COMS DONE')
  # sys.exit(0)
  
  # Want to write unwrapped coordinates?
  if (qcom is True):
    print(f'# Write unwrapped COM of individual molecules')
    write_com_mol(traj_com_unwrap,nmol,framenum,odir)
  
  print(f'# MSD calculations')
  nfreq = int(framenum/interval)
  taus = []; [taus.append(interval*i) for i in range(0,nfreq)];
  ntau = len(taus) # number of data points along the delay time
  
  # Setup msd arrays for individual molecule types
  msd = mymsd.setup_msd_arrays(ntype,ntau)
  
  # Calculate MSD for delay times, tau in {taus}
  mymsd.calculate_msd(taus,framenum,traj_com_unwrap,id_type,msd)
  
  ### SP: NEED to make something consistent with STANALYSER output path ??!!
  # Write MSD outputs
  write_msd_outputs(msd,taus,ntype,name_type,odir)
  
  return



def main(settings: Optional[dict] = None) -> None:

    if settings is None:
         settings = dict(sta.get_settings(ANALYSIS_NAME))

    # non-system arguments will be handled at the beginnig of this function
    run_msd_solution(**settings)



if __name__ == '__main__':

    main()



