import argparse
import sys
from typing import Optional, cast

#import stanalyzer.bin.stanalyzer as sta
import MDAnalysis as mda    # type: ignore
import numpy as np

from MDAnalysis.analysis import distances # type: ignore

ANALYSIS_NAME = 'contacts'

def header(outfile: Optional[sta.FileLike] = None, np_formatted=False) -> str:
    """Returns a header string and, if optionally writes it to a file

    If np_formatted is true, the `#` is omitted."""
    if np_formatted:
        header_str = "time contacts"
    else:
        header_str = "#time contacts"

    print(header_str, file=outfile)

    return header_str
    
def write_contacts(psf: sta.FileRef, traj: sta.FileRefList, sel: str,
               out: sta.FileRef, contact_threshold: float = 5.0,ref_psf: Optional[sta.FileRef] = None,
               ref_coor: Optional[sta.FileRef] = None,
               ref_frame_type: str = 'specific', 
               ref_frame_num: int = 1,
               interval: int = 1) -> None:
    """Writes contacts to `out` file"""

# Parameters
topology_file = 'topology.pdb'   # Path to the topology file
trajectory_file = 'trajectory.dcd'  # Path to the trajectory file
contact_threshold = 5.0  # Contact threshold distance in Angstroms

# Load the universe
u = mda.Universe('topology.pdb', 'trajectory.dcd')

# Function to calculate contacts
def calculate_contacts(universe, threshold):
    # Select all residues
    residues = u.select_atoms('all')
    
    # Create a list to store contact information
    contact_residues = set()
    
    # Iterate over frames in the trajectory
    for ts in universe.trajectory:
        # Select all atoms
        atoms = u.select_atoms('all')
        
        # Calculate the pairwise distances
        pairwise_distances = distances.distance_array(
            atoms.positions, atoms.positions, box=atoms.dimensions
        )
        
        # Check distances between residues
        num_residues = len(residues)
        for i in range(num_residues):
            for j in range(i + 1, num_residues):
                # Get the residue positions
                residue_i = residues[i].atoms
                residue_j = residues[j].atoms
                
                # Calculate distance between the center of mass of residues
                com_i = residue_i.center_of_mass()
                com_j = residue_j.center_of_mass()
                
                # Calculate distance
                dist = np.linalg.norm(com_i - com_j)
                
                if dist < threshold:
                    contact_residues.add((residues[i].resname, residues[i].resid, residues[j].resname, residues[j].resid))
    
    return contact_residues
# Calculate contacts
contacts = calculate_contacts(u, contact_threshold)

def get_parser() -> argparse.ArgumentParser:
   parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
   sta.add_project_args(parser, 'pdb', 'traj', 'contact_threshold', 'out', 'interval')
   parser.add_argument('--contact_threshold', type=float, metavar='N', default='5.0',
                       help="Contact threshold distance in Angstroms.") 
   parser.add_argument('--sel', metavar='selection',
                       help="Atom selection for contacts calculation")
   parser.add_argument('-c', '--center', action='store_true')
   parser.add_argument('--sel', metavar='selection')

   return parser

def main():
   u = mda.Universe('topology.pdb', 'trajectory.dcd')
 
def main(settings: Optional[dict] = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_contacts(**settings)


if __name__ == '__main__':
    main()

