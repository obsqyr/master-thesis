#!/usr/bin/env python3
from ase.io.vasp import *
import numpy as np
import os

# Global variables
infiles = ['forces.infile', 'positions.infile', 'potentials.infile', 'numatoms.infile']
    
def vasp_read(directory, it):
    '''
    Reads from a VASP directory and produces a data set of
    positions, forces and potentials (eventually more).

    Parameters:
    directory (str): name of the targeted directory
    it (int): number of timesteps to be included
    '''
    print("reading VASP")
    
    sl = slice(0, None)
    atoms = list(read_vasp_xml(filename=directory + 'vasprun.xml', index=sl))
    first = True
    for i, atom in enumerate(atoms):
        if first:
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/', True)
            first = False
        else:
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/')
            
            
def write_atom_to_infiles(atoms, directory, num=False):
    '''
    Writes an atom object to files. Writes forces, positions, 
    potential energy for each timestep and finally number of atoms.
    '''
    if not os.path.isdir(directory):
        os.mkdir(directory)
        
    if not num:
        with open(directory+'forces.infile', 'a+') as f:
            np.savetxt(f, atoms.get_forces(), fmt='%-1.7f')
        with open(directory+'positions.infile', 'a+') as f:
            np.savetxt(f, atoms.get_positions(), fmt='%-1.7f')
        with open(directory+'potentials.infile', 'a+') as f:
            f.write(str(atoms.get_potential_energy()) + '\n')
    else:
        with open(directory+'numatoms.infile', 'a+') as f:
            f.write(str(len(atoms.get_chemical_symbols()))+'\n')
            f.write(str(atoms.get_atomic_numbers()[0]))
            

def clear_infiles(directory):
    '''
    Delete current infiles in directory
    '''
    #infiles = ['forces.infile', 'positions.infile', 'potentials.infile']
    for f in infiles:
        if os.path.exists(directory[:-1]+'_infiles/'+f):
            print("deleting " + directory[:-1]+'_infiles/'+f)
            os.remove(directory[:-1]+'_infiles/'+f)
    
def read_infiles(directory):
    dir = directory[:-1]+'_infiles/'
    
    forces = np.loadtxt(dir+'forces.infile')
    positions = np.loadtxt(dir+'positions.infile')
    potentials = np.loadtxt(dir+'potentials.infile')
    num_atoms = np.loadtxt(dir+'numatoms.infile', dtype=int)
    
    return forces, positions, np.array(potentials), num_atoms

if __name__ == "__main__":
    clear_infiles("Al_300K/")
    vasp_read("Al_300K/")
    #f, pos, pot, num = read_infiles("Al_300K/")
