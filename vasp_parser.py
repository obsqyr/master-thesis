#!/usr/bin/env python3
from ase.io.vasp import *
import numpy as np
import os
from ase.io.trajectory import Trajectory
import properties as pr
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

# Global variables
infiles = ['forces.infile', 'positions.infile', 'potentials.infile', 'numatoms.infile', 'velocities.npy']
    
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
        
    # delete trajectory file if it exits
    if os.path.exists(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj'):
        os.remove(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj')

    # write atoms to trajectory file
    traj = Trajectory(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj',mode='a', atoms=atoms)

    # write relevant data to file
    first = True
    velocities = []
    for i, atom in enumerate(atoms):
        if first:
            prev_atom = atom
            traj.write(atom)
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/', True)
            first = False
        else:
            # only support for cubic unit cells
            unit_cell = atom.get_cell()
            print('unit_cell', unit_cell)
            #print('masses', atom.get_masses())
            #print('unit_cell', list(unit_cell)[0][0])
            l = list(unit_cell)[0][0]
            if l == 0:
                l = list(unit_cell)[0][1]
            print('l', l)

            # convert positions to fractional form
            pos = atom.get_positions()
            old_pos = prev_atom.get_positions()
            pos = np.array([i/l for i in pos])
            old_pos = np.array([i/l for i in old_pos])
            
            #v = atom.get_positions() - prev_atom.get_positions()
            v = pos - old_pos # unit: 1 / fs
            v = [pr.normalize_half(i,l) for i in v]
            #print('v normalized', np.array(v))
            #print('atom pos', atom.get_positions())
            #print('prev_atom pos', prev_atom.get_positions())
            # does this factor of 10 make sense?
            v = np.array([i*l for i in v]) # unit: Å / fs
            print('v in Å', v)
            velocities.append(v)
            prev_atom = atom
            
            traj.write(atom)
            write_atom_to_infiles(atom, directory[:-1]+'_infiles/')
    print(len(velocities))
    velocities = np.array(velocities)
    #print(velocities[:10])
    # write velocities to file
    np.save(directory[:-1]+'_infiles/'+'velocities.npy', velocities)
        
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

def calculate_properties_vasp(element, eq):
    '''
    Calculates properties of atoms trajectory objects. Writes to
    properties file. Requires that the trajectory file has 
    been read from vasp. 
    '''
    # generate trajectory file, reading from vasp data
    
    # import atoms from trajectory file
    traj = Trajectory(element +'_300K_infiles/'+ element +'.traj')
    atoms = [atom for atom in traj]
    id = str(element)+'_'+'DFT_eq_' + str(eq)
    
    # delete old DFT properties file, it it exists
    if os.path.exists('property_calculations/properties_' + id + '.txt'):
        os.remove('property_calculations/properties_' + id + '.txt')

    pr.initialize_properties_file(atoms[1], atoms[0][:4], id, 5, True) 
    #for atom in atoms:
    print('atoms', atoms[0])
    velocities = np.load(element +'_300K_infiles/velocities.npy')
    #print(velocities[:10])
    #print(velocities[0])
    #for atom in atoms:
    #    print('before', atom.get_velocities())
    #    MaxwellBoltzmannDistribution(atom, temperature_K = 300)
    #    print('after', atom.get_velocities())
    MSDs = []
    Cvs = []
    temps = []
    
    for i in range(len(atoms) - eq):
        i += eq
        # add velocities for all atoms but the last (since the
        # last atom is not going to be moved)
        if i != len(atoms) - 1:
        #    print(i)
            atoms[i].set_velocities(velocities[i])
            #print(atoms[i].get_velocities())
        #print(atom.get_forces())
        #print(atoms[i].get_positions())
        #print(atoms[i].get_momenta())
        if i % 100 == 0:
            #print(i)
            _, _, _, t = pr.energies_and_temp(atoms[i])
            temps.append(t)
            Cvs.append(pr.specific_heat(temps, len(atoms[i]), atoms[i]))
            MSDs.append(pr.meansquaredisp(atoms[i], atoms[eq]))
            pr.calc_properties(atoms[eq], atoms[i], id, 5, True)
    pr.finalize_properties_file(atoms[-1], id, 5, True, True)

    # show whether or not MSD converges
    MSD_averaged = []
    MSD = 0
    for i, m in enumerate(MSDs):
        if i != 0:
            MSD += m
            MSD_averaged.append(MSD / i)
            #print('MSD / i', MSD / i)
            #print(MSD/(i), MSD, i)
    print('MSD_averaged', len(MSD_averaged))
    
    # show whether or not specific heat converges
    Cvs_averaged = []
    Cv = 0
    for i, c in enumerate(Cvs):
        Cv += c
        Cvs_averaged.append(Cv / (i+1))

    return MSD_averaged, Cvs_averaged
    
if __name__ == "__main__":
    clear_infiles("Al_300K/")
    #clear_infiles("Si_300K/")
    
    vasp_read("Al_300K/", "xml")
    #vasp_read("Si_300K/", "OUTCAR")
    #f, pos, pot, num = read_infiles("Al_300K/")
    #read_vasp_out("Si_300K/OUTCAR")    
    #calculate_properties_vasp('Si', 0)
    #calculate_properties_vasp('Al', 0)
    #calculate_properties_vasp('Si', 6000)
    calculate_properties_vasp('Al', 0)
    #calculate_properties_vasp('Al', 2000)
