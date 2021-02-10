#!/usr/bin/env python3
from ase.io.vasp import *
import numpy as np
import os

# Global variables
infiles = ['forces.infile', 'positions.infile', 'potentials.infile']
    
def vasp_read(directory, it):
    '''
    Reads from a VASP directory and produces a data set of
    positions, forces and potentials (eventually more).

    Parameters:
    directory (str): name of the targeted directory
    it (int): number of timesteps to be included
    '''

    print("reading VASP")
    #print(read_vasp(directory + "POSCAR"))
    #print(directory + "OUTCAR")
    #read_vasp_out(directory + "OUTCAR", index=0)

    for i in range(0, it):
        print(i)
        atoms = [d for d in read_vasp_xml(filename=directory + 'vasprun.xml', index=i)][0]
        write_atoms_to_infile(atoms)
        
    #print(vasp_data)
    #print(vasp_data.get_potential_energy())
    
    #print(read_vasp_xdatcar(directory + "XDATCAR")[0])

def write_atoms_to_infile(atoms):
    '''
    Writes an atom object to files. Writes forces, positions and 
    potential energy.
    '''
    with open('forces.infile', 'a+') as f:
        np.savetxt(f, atoms.get_forces(), fmt='%-1.7f')
    with open('positions.infile', 'a+') as f:
        np.savetxt(f, atoms.get_positions(), fmt='%-1.7f')
    with open('potentials.infile', 'a+') as f:
        f.write(str(atoms.get_potential_energy()) + '\n')
        #np.savetxt(f, atoms.get_potential_energy(), fmt='%-1.7f')

def clear_infiles():
    '''
    Delete current infiles.
    '''
    #infiles = ['forces.infile', 'positions.infile', 'potentials.infile']
    for f in infiles:
        if os.path.exists(f):
            print("deleting " + f)
            os.remove(f)
    
def read_infiles():
    forces = np.loadtxt('forces.infile')
    positions = np.loadtxt('positions.infile')
    
    potentials = []
    with open('potentials.infile', 'r') as f:
        for line in f:
            potentials.append(float(line))

    return forces, positions, np.array(potentials)

if __name__ == "__main__":
    clear_infiles()
    vasp_read("Al_300K/", 100)
    #f, pos, pot = read_infiles()
    
