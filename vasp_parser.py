#!/usr/bin/env python3
from ase.io.vasp import *
import numpy as np
import os
from ase.io.trajectory import Trajectory
import properties as pr
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


# Global variables
infiles = ['forces.infile', 'positions.infile', 'potentials.infile', 'numatoms.infile']
    
def vasp_read(directory, filetype):
    '''
    Reads xml or OUTCAR file from a VASP directory and 
    produces a data set of positions, forces and 
    potentials (eventually more) and a trajectory file.

    Parameters:
    directory (str): name of the targeted directory
    filetype (str): whether xml or OUTCAR file should be read
    '''
    print("reading VASP "+ filetype +" from " + directory)
    
    sl = slice(0, None)
    # read correct filetype
    if filetype == "xml":
        atoms = list(read_vasp_xml(filename=directory + 'vasprun.xml', index=sl))
    elif filetype == "OUTCAR":
        atoms = list(read_vasp_out(directory + 'OUTCAR', index=sl))
    else:
        raise TypeError('Invalid filetype to read')
        
    # write atoms to trajectory file
    traj = Trajectory(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj',mode='a', atoms=atoms)

    # write relevant data to file
    first = True
    for i, atom in enumerate(atoms):
        if first:
            traj.write(atom)
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/', True)
            first = False
        else:
            traj.write(atom)
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/')
            
def write_atom_to_infiles(atoms, directory, num=False):
    '''
    Writes an atom object to files. Writes forces, positions, 
    potential energy for each timestep and finally number of atoms.
    '''
    if not os.path.isdir(directory):
        os.mkdir(directory)

    with open(directory+'forces.infile', 'a+') as f:
        np.savetxt(f, atoms.get_forces(), fmt='%-1.7f')
    with open(directory+'positions.infile', 'a+') as f:
        np.savetxt(f, atoms.get_positions(), fmt='%-1.7f')
    with open(directory+'potentials.infile', 'a+') as f:
        f.write(str(atoms.get_potential_energy()) + '\n')
    if num:
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


def calculate_properties_vasp(element):
    '''
    Calculates properties of atoms trajectory objects. Requires
    that the trajectory file has been read from vasp. 
    '''
    # generate trajectory file, reading from vasp data
    
    # import atoms from trajectory file
    traj = Trajectory(element +'_300K_infiles/'+ element +'.traj')
    atoms = [atom for atom in traj]
    id = str(element)+'_'+'DFT'
    
    # delete old DFT properties file, it it exists
    if os.path.exists('property_calculations/properties_' + id + '.txt'):
        os.remove('property_calculations/properties_' + id + '.txt')

    pr.initialize_properties_file(atoms[1], atoms[0][:4], id, 5, True) 
    #for atom in atoms:
    print('atoms', atoms[0])
    for atom in atoms:
        MaxwellBoltzmannDistribution(atom, temperature_K = 300)
    for i, atom in enumerate(atoms):
        #print(atom.get_forces())
        if i % 100 == 0:
            print(i)
            pr.calc_properties(atoms[i-1], atoms[i], id, 5, True)
    pr.finalize_properties_file(atoms[-1], id, 5, True, True)

if __name__ == "__main__":
    #clear_infiles("Al_300K/")
    #clear_infiles("Si_300K/")
    
    #vasp_read("Al_300K/", "xml")
    #vasp_read("Si_300K/", "OUTCAR")
    #f, pos, pot, num = read_infiles("Al_300K/")
    #read_vasp_out("Si_300K/OUTCAR")    
    calculate_properties_vasp('Al')
