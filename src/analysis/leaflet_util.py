#!/usr/bin/python
import sys
import MDAnalysis as mda
import MDAnalysis.analysis.leaflet as mdaleaflet
import numpy as np

nside=2 # number of leaflets in a bilayer !!!

def find_min_dist(pos,pos_grp):
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

  ngrp = len(pos_grp) # number of vectors in the position group

  min_dist = np.linalg.norm(pos-pos_grp[0])  # initial dist
  if (ngrp > 1):
    for i in range(1,ngrp):
      tdist = np.linalg.norm(pos-pos_grp[i,:])
      if (tdist < min_dist): min_dist = tdist

  return (min_dist)



def assign_groups_to_leaflets(sel_grp,ngroups):
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
  ng=np.zeros([ngroups],dtype=int)
  pos_grp=[] # position vector for sel_grp
  for i in range(0,ngroups):
    tmp_grp = sel_grp[i]
    ng[i] = len(tmp_grp)
    pos_grp.append(tmp_grp.positions)

  # strat from the 3rd group
  for i in range(2,ngroups):
    if (ng[i] < int(0.1*ng[1]) and ng[i] < int(0.1*ng[0])):
      # Merge members into a closer group between the largest two
      for j in range(0,ng[i]):
        # find min. distances to group0 and group 1
        tpos=pos_grp[i][j,:]
        tdmin0 = find_min_dist(tpos,pos_grp[0])
        tdmin1 = find_min_dist(tpos,pos_grp[1])

        # now compare two min dist
        if (tdmin0 < tdmin1): # put this to B.groups(0)
          # print(f'assignment {j} to group0')
          sel_grp[0]+=sel_grp[i].atoms[j]
        else: # put this to B.groups(1)
          # print(f'assignment {j} to group1')
          sel_grp[1]+=sel_grp[i].atoms[j]
      # print("# assignment done")
  return (sel_grp)



def assign_leaflet_zpos(u,atomgroup):
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

  leaflets = []
  pos = atomgroup.positions # position vector
  zmin = np.min(pos[:,2])
  zmax = np.max(pos[:,2])
  # set z value to separate leaflet
  zcent = 0.5 * (zmin + zmax)

  # upper leaflet; z > zcent
  tag0 = atomgroup.intersection(u.select_atoms("prop z > %g" % zcent))
  leaflets.append(tag0)

  # lower leaflet; z < zcent
  tag1 = atomgroup.intersection(u.select_atoms("prop z < %g" % zcent))
  leaflets.append(tag1)

  # remainder - neither upper nor lower leaflet
  tag = atomgroup.difference(tag0+tag1)
  leaflets.append(tag)
  
  return (leaflets)



def assign_leaflet_mda(u,atomgroup):
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

  B = mdaleaflet.LeafletFinder(u,atomgroup)
  leaflets=[]
  for i in range(0,len(B.groups())):
    leaflets.append(B.groups(i))
    # print(f'simplest: len(B.groups({i}))= {len(leaflets[i])}')
    # if(len(leaflets[i]) < 10): print(leaflets[i])

  return (leaflets)



def assign_leaflet_simple(u,atomgroup):
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

  dcut=15.0      # default from MDAnalysis
  dmin,dmax=(0.0,2.0*dcut)

  # Start with default setting
  c_dcut=0.5*(dmin+dmax)
  B = mdaleaflet.LeafletFinder(u,atomgroup,cutoff=c_dcut,pbc=True)
  ngroups=len(B.groups())

  # a hybrid stride-bisection method to obtain two groups
  if (ngroups != 2):
    # initial direction
    if (ngroups > 2):        # increase dmax by the width of [dmin,dmax]
      dmax += (dmax-dmin)
    else:                    # decrease dmax by the half-widht of [dmin,dmax]
      dmax -= 0.5*(dmax-dmin)

    ## print(f'# iter 0, len(B.groups)={len(B.groups()} dcut={c_dcut}')

    imax=100 # max iteration
    for i in range(1,imax+1):
      # save len(B.groups()) and cutoff distance from the previous iteration
      p_ngroups=ngroups
      p_dcut=c_dcut

      # set the current cutoff distance
      c_dcut=0.5*(dmin+dmax)
      # assign leaflet with c_dcut
      B = LeafletFinder(u,atomgroup,cutoff=c_dcut,pbc=True)
      ngroups=len(B.groups())

      # check the current len(B.groups())
      if (ngroups == 2):
        ## print(f'# simple iter {i}, len(B.groups)={ngroups} after LeafletFinder, dcut= {c_dcut}')
        ## print("# B.groups", B.groups())
        break
      else:
        if (ngroups < 2): # decrease the cutoff distance for the next iteration
          if (p_ngroups > 2): # p_dcut is too small
            if (p_dcut > c_dcut):
               print(f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & c_dcut {c_dcut} for cgr= {ngroups}')
               sys.exit(1)
            # update dmin,dmax 
            dmin,dmax = (p_dcut,c_dcut)
          else: # move along the same direction, decrease dmax
            dmax -= 0.5*(dmax-dmin) # decrease dmax by the half-width of [dmin,dmax]
        elif (ngroups > 2): # increase the cutoff distance for the next iteration
          if (p_ngroups < 2): # p_dcut is too large
            if (p_dcut < c_dcut):
               print(f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & c_dcut {c_dcut} for cgr= {ngroups}')
               sys.exit(1)
            # update dmin,dmax 
            dmin,dmax = (c_dcut,p_dcut)
          else: # move along the same direction, increase dmax
            dmax += (dmax-dmin) # increase dmax by the width of [dmin,dmax]

        ## print(f'# simple: iter {i}, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

    if (ngroups != nside):
      print(f'# Job failed! Exit script !!')
      for j in range(0,len(B.groups())):
        if(len(B.groups(j))<10): print(B.groups(j))
      sys.exit(1)

  # return assigned atom groups
  leaflets=[]
  for i in range(0,nside):
    leaflets.append(B.groups(i))
    ## print(f'# simple: Lflt{i}: len(B.groups({i})) = {len(leaflets[i])}')

  return (leaflets)



def assign_leaflet(u,atomgroup):
  """
  ----------
  Assign leaflet using a hybdird stride-bisection method when necessary.
  ----------
  NOTE: IN typical situations that chosen atoms are appropriate for lipid tails,
        leaflets should have comparable numbers of atoms.

  input
        u        : MDA Universe
        atomgroup: atom group to be assinged

  output
        leaflets : atom groups for individual leaflets
  """


  dcut=15.0      # default from MDAnalysis
  dmin,dmax=(0.0,2.0*dcut)

  # Start with default setting
  c_dcut=0.5*(dmin+dmax)
  B = mdaleaflet.LeafletFinder(u,atomgroup,cutoff=c_dcut,pbc=True)
  ngroups=len(B.groups())
  # 
  sel_grp,ng_sel = [],[]
  for ig in range(0,ngroups):
    sel_grp.append(B.groups(ig))
    ng_sel.append(len(sel_grp[ig]))

  # check if selected groups for the leaflets are reasonable
  if (ngroups == 2):
    if(ng_sel[1] > int(0.5*ng_sel[0])): qdone = True
    else: qdone = False # suspicious case; more groups may be assigned with smaller dcut
  elif (ngroups > 2):
     # special case: the 3rd group has significantly smaller members than the 2nd group
     if (ng_sel[2] < int(0.1*ng_sel[1]) and \
         ng_sel[2] < int(0.1*ng_sel[0])):
       sel_grp = assign_groups_to_leaflets(sel_grp,ngroups)
       qdone = True
     else: qdone = False
  else: qdone = False

  # a hybrid stride-bisection method to obtain two groups
  ## if (ngroups != 2 and not qdone):
  if (not qdone):
    # initial direction
    if (ngroups == 2):       # suspicious case that needs decreased cutoff dist.
      dmax -= 0.5*(dmax-dmin)
    elif (ngroups > 2):      # increase dmax by the width of [dmin,dmax]
      dmax += (dmax-dmin)
    else:                    # decrease dmax by the half-widht of [dmin,dmax]
      dmax -= 0.5*(dmax-dmin)

    ## print(f'# iter 0, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

    imax=100 # max iteration
    for i in range(1,imax+1):
      # save len(B.groups()) and cutoff distance from the previous iteration
      p_ngroups=ngroups
      p_dcut=c_dcut

      # set the current cutoff distance
      c_dcut=0.5*(dmin+dmax)
      # assign leaflet with c_dcut
      B = mdaleaflet.LeafletFinder(u,atomgroup,cutoff=c_dcut,pbc=True)
      ngroups=len(B.groups())
      sel_grp,ng_sel = [],[]
      for ig in range(0,ngroups):
        sel_grp.append(B.groups(ig))
        ng_sel.append(len(sel_grp[ig]))

      # check the current iteration
      if (ngroups == 2):
        if( len(sel_grp[1]) > int(0.5*len(sel_grp[0]))): qdone = True
        else: qdone = False

        if (not qdone):
          ## suspicious case, need to decrease the cutoff distance
          ## i.e., other lipids vs. chol in the middle of flip-floping 
          dmax -= 0.5*(dmax-dmin) # decrease dmax by the half-width of [dmin,dmax]
        else:
          ## print(f'# iter {i}, len(B.groups) = {ngroups} after LeafletFinder, dcut= {c_dcut}')
          ## print(f'# B.groups', B.groups())
          break
      else:
        if (ngroups > 2): # increase the cutoff distance for the next iteration
          # special case:
          # the 3rd group has significantly smaller members than the 2nd group
          if (ng_sel[2] < int(0.1*ng_sel[1]) and \
               ng_sel[2] < int(0.1*ng_sel[0])):
            sel_grp = assign_groups_to_leaflets(sel_grp,ngroups)
            qdone = True
          else: qdone = False
          if (qdone): break

          if (p_ngroups < 2): # p_dcut is too large
            if (p_dcut < c_dcut):
               print(f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & c_dcut {c_dcut} for cgr= {ngroups}')
               sys.exit(1)
            # update dmin,dmax 
            dmin,dmax = (c_dcut,p_dcut)
          else: # move along the same direction, increase dmax
            dmax += (dmax-dmin) # increase dmax by the width of [dmin,dmax]

        elif (ngroups < 2): # decrease the cutoff distance for the next iteration
          if (p_ngroups > 2): # p_dcut is too small
            if (p_dcut > c_dcut):
               print(f'# errorneous p_dcut {p_dcut} for pgr= {p_ngroups} & c_dcut {c_dcut} for cgr= {ngroups}')
               sys.exit(1)
            # update dmin,dmax 
            dmin,dmax = (p_dcut,c_dcut)
          else: # move along the same direction, decrease dmax
            dmax -= 0.5*(dmax-dmin) # decrease dmax by the half-width of [dmin,dmax]

        ## print(f'# iter {i}, len(B.groups) = {len(B.groups())} dcut={c_dcut}')

    # if (ngroups != nside):
    if (qdone): pass
    else:
      print(f'# Job failed! Exit script !!')
      print(len(B.groups()))
      for j in range(0,len(B.groups())):
        if (len(B.groups(j))<10): print(B.groups(j))
        # print(B.groups(j))
      sys.exit(1)

  ## print(f'# assign_leaflet: len(B.groups()) = {ngroups}')

  # return assigned atom groups
  leaflets=[]
  for i in range(0,nside):
    leaflets.append(sel_grp[i])
    ## print(f'# assgin_leaflet: leaflet {i}: nlipids = {len(leaflets[i])}')

  return (leaflets)


