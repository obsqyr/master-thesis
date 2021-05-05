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
        
    # delete trajectory file if it exits
    if os.path.exists(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj'):
        os.remove(directory[:-1]+'_infiles/'+atoms[0].get_chemical_symbols()[0]+'.traj')

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
    for atom in atoms:
        MaxwellBoltzmannDistribution(atom, temperature_K = 300)
    MSDs = []
    Cvs = []
    temps = []
    
    ### TEST MSD
    x0 = [[0,      0,      0     ],
          [0,      2.02395, 2.02395],
          [2.02395, 0,      2.02395],
          [2.02395, 2.02395, 0     ],
          [0,      0,      4.0479 ],
          [0,      2.02395, 6.07185],
          [2.02395, 0,      6.07185],
          [2.02395, 2.02395, 4.0479 ],
          [0,      4.0479,  0     ],
          [0,      6.07185, 2.02395],
          [2.02395, 4.0479,  2.02395],
          [2.02395, 6.07185, 0     ],
          [0,      4.0479,  4.0479 ],
          [0,      6.07185, 6.07185],
          [2.02395, 4.0479,  6.07185],
          [2.02395, 6.07185, 4.0479 ],
          [4.0479,  0,      0     ],
          [4.0479,  2.02395, 2.02395],
          [6.07185, 0,      2.02395],
          [6.07185, 2.02395, 0     ],
          [4.0479,  0,      4.0479 ],
          [4.0479,  2.02395, 6.07185],
          [6.07185, 0,      6.07185],
          [6.07185, 2.02395, 4.0479 ],
          [4.0479,  4.0479,  0     ],
          [4.0479,  6.07185, 2.02395],
          [6.07185, 4.0479,  2.02395],
          [6.07185, 6.07185, 0     ],
          [4.0479,  4.0479,  4.0479 ],
          [4.0479,  6.07185, 6.07185],
          [6.07185, 4.0479,  6.07185]]
    x0 = np.array(x0)
    
    x_last = [[ 0.0422451,   0.01912541, -0.05709584],
              [-0.02571535,  1.89871802,  2.00477168],
              [ 2.11881617,  0.06492702,  1.9534501 ],
              [ 2.13001194,  2.05870131, -0.0852118 ],
              [-0.06088588, -0.08552464,  3.9965635 ],
              [-0.14397687,  1.88148807,  6.01665408],
              [ 2.0066865,   0.01318263,  6.03086966],
              [ 2.02498751,  1.94527234,  4.07938469],
              [ 0.07118548,  3.92718918, -0.07103578],
              [ 0.07878729,  6.21525013,  2.14018978],
              [ 2.04456109,  4.02217587,  2.10467059],
              [ 2.07168461,  5.95612501, -0.01349559],
              [ 0.08744669,  4.02104453,  3.9405703 ],
              [ 0.10457843,  6.04533795,  6.1188013 ],
              [ 1.98416628,  4.02548036,  6.00015813],
              [ 2.00876304,  6.12576043,  4.04940943],
              [ 3.91748002,  0.03355322, -0.10687896],
              [ 4.02287771,  2.00819978,  2.04070709],
              [ 5.97843322,  0.09659275,  2.07977856],
              [ 5.96922893,  2.10123714, -0.07469658],
              [ 4.07382142, -0.11106306,  4.19304867],
              [ 4.12156885,  2.15604782,  6.0372239 ],
              [ 6.04683152, -0.06115595,  6.19747028],
              [ 6.0484283,   2.10079206,  4.05411977],
              [ 4.09074193,  4.135346,    0.15642483],
              [ 4.11440776,  6.11104686,  2.05719678],
              [ 6.07859526,  3.92643561,  2.02595015],
              [ 6.06041137,  6.17509882,  0.13202432],
              [ 4.06001182,  4.19188633,  4.03420771],
              [ 3.94435263,  6.07252489,  6.13097621],
              [ 6.06726869,  4.13728889,  6.0240104 ],
              [ 6.02736259,  6.04432941,  3.94577905]]
    x_last = np.array(x_last)

    x0 = atoms[2000].get_positions()
    x_last = atoms[9000].get_positions()

    def distance2(pos1, pos2):
        """Calculates the sqared distance between two atomsx in 3D space"""
        return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2

    #print(x0, x_last)
    msd = 0
    for i in range(len(x0)): 
        #print('distance^2 np',np.linalg.norm(x0[i] - x_last[i])**2)
        #print('distance^2 d2', distance2(x0[i], x_last[i]))
        msd += np.linalg.norm(x0[i] - x_last[i])**2
        
    #print('MSD', msd/len(x0))

    for i in range(len(atoms) - eq):
        i += eq
        #print(atom.get_forces())
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
    clear_infiles("Si_300K/")
    
    vasp_read("Al_300K/", "xml")
    vasp_read("Si_300K/", "OUTCAR")
    #f, pos, pot, num = read_infiles("Al_300K/")
    #read_vasp_out("Si_300K/OUTCAR")    
    calculate_properties_vasp('Si', 0)
    calculate_properties_vasp('Al', 0)
    calculate_properties_vasp('Si', 2000)
    calculate_properties_vasp('Al', 2000)
