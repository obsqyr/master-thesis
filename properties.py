#!/usr/bin/env python3

import math
from ase import Atoms
from ase import units
import numpy as np
import os
import chemparse
from ase.build import bulk
import sklearn

from read_settings import read_settings_file

# This file contains functions to calculate material properties

def normalize_half(vec):
    '''
    input:
    np.array of dimension 3 (position vector)
    Add/remove an integer +/-N to each element to place it in 
    the range [-1/2,1/2)
    This is useful to find the shortest vector C between two 
    points A, B in a space with periodic boundary conditions [0,1):
           C = (A-B).normalize_half()
    '''
    '''
    temp = []
    for v in vec:
        if v > 0.5:
            temp.append(v-1)
        elif v < -0.5:
            temp.append(v+1)
        else:
            temp.append(v)
    return np.array(temp)
    '''
    
    #vec_std = (vec - np.amin(vec)) / (np.amax(vec) - np.amin(vec))
    
    #vec_std = (vec - np.amin(vec)) / 
    #vec_scaled = vec_std * (1/2 - (-1/2)) + -1/2

    #scale = (1/2 - (-1/2)) / (np.amax(vec) - np.amin(vec))
    #vec_scaled2 = scale * vec + np.amin(vec) - np.amin(vec) * scale
    #print(vec_scaled)
    #print(vec_scaled2)
    #x - (x*2+1)//2

    return vec - (vec*2+1)//2 

def specific_heat(temp_store, N, atoms):
    """Calculates the specific heat for a material.
    Given by the formula: (E[T²] - E[T]²)/ E[T]² = 3*2^-1*N^-1*(1-3*N*Kb*2^-1*Cv^-1).
    Where Kb is boltzmansconstant, N is the total number of atoms, T is temperature and Cv the specific heat.
    E[A(t)] calculates the expectation value of A, which can in this case be seen as a time average for the
    phase variable A(t).
    Parameters:
    temp_store (list): The list over all intantaneous temperatures of a material once MD has been run.
    N (int): The total number of atoms in the material.
    Returns:
    float: specific heat is returned (J/(K*Kg))
    """
    if len(temp_store) == 0:
        raise ValueError("temp_store is empty, invalid value.")
    steps = len(temp_store)
    z = sum(atoms.get_masses()) * units._amu # total mass: atomic units to kg
    # Set M = (E[T²] - E[T]²)/ E[T]²
    ET = sum(temp_store)/steps
    ET2 = sum(np.array(temp_store)**2)/steps
    M = (ET2 - ET**2)/ET**2
    settings = read_settings_file()
    N = N / settings['supercell_size']**3
    Cv = ((9*ET**2*N*units._k) / (ET**2 * (6+4*N) - 4*N*ET2)) / z * settings['supercell_size']**3
    return Cv

def specific_heat_NVT(energies, N, atoms, T=300):
    """Calculates the specific heat for a material.
    energies (list): The list over all intantaneous total energies of a material once MD has been run.
    N (int): The total number of atoms in the material.
    Returns:
    float: specific heat is returned (J/(K*Kg))
    """
    if len(energies) == 0:
        raise ValueError("energies is empty, invalid value.")
    steps = len(energies)
    z = sum(atoms.get_masses()) * units._amu # total mass: atomic units to kg
    # convert energies from eV to J
    energies = [i*units._e*N for i in energies]
    
    EE = sum(energies)/steps
    EE2 = sum(np.array(energies)**2)/steps
    settings = read_settings_file()
    #Cv = ((9*ET**2*N*units._k) / (ET**2 * (6+4*N) - 4*N*ET2)) / z * settings['supercell_size']**3
    #T = 300 # assume T always 300 K, try instantaneous temperature also
    settings = read_settings_file()
    Cv = (EE2 - EE**2) / (units._k * T**2) / z
    return Cv

def distance2(pos1, pos2):
    """Calculates the sqared distance between two atomsx in 3D space"""
    return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2

def distance(pos1, pos2):
    return math.sqrt(distance2(pos1, pos2))

def meansquaredisp(atoms, old_atoms):
    """ Calculates the mean squared displacement
    Parameters:
    atoms (obj):atoms is an atom object from ase.
    old_atoms (obj):old_atoms is an atom object from the python library.
    Returns:
    float: The mean squared displacement.
    """
    pos = atoms.get_positions()
    old_pos = old_atoms.get_positions()
    length = len(pos)
    
    # only support for cubic unit cells
    unit_cell = atoms.get_cell()
    #print('unit_cell', list(unit_cell)[0][0])
    l = list(unit_cell)[0][0]
    if l == 0:
        l = list(unit_cell)[0][1]
    #print('l', l)

    # convert positions to fractional form
    pos = np.array([i/l for i in pos])
    old_pos = np.array([i/l for i in old_pos])
    #print('pos', pos)
    #print('old_pos', old_pos)
    
    #np.set_printoptions(suppress=True)
    #print('old_pos', old_pos, 'pos', pos)
    if length != len(old_pos):
        raise TypeError("Number of atoms doesnt match.")
        sys.exit('ERROR')

    msd = 0.0
    #msd_np = 0.0
    for atom in range(length):
        #msd += distance2(pos[atom], old_pos[atom])
        diff = pos[atom] - old_pos[atom]
        #print(diff)
        # do I need to do this?
        #diff = [i/l for i in diff]
        diff_norm = normalize_half(diff)
        #print('diff_norm', diff_norm)
        msd += np.linalg.norm(diff_norm)**2
        #print(msd)

    #print('MSD', msd/length)
    return msd/length

def energies_and_temp(a):
    """ Calculates the energies and temperature.
    Parameters:
    a (obj): a is an atoms object of class defined in ase.
    Returns:
    tuple: returns a tuple of potential energi, kinetic energy, total energy
            and temperature.
    """
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    etot = epot + ekin
    t = ekin / (1.5 * units.kB)
    return epot, ekin, etot, t

def lattice_constants(a):
    """ Calculates the lattice constant of a materialself.
    Parameters:
    a (obj): a is an atoms object of class defined in ase.
    Returns:
    list: returns the lattice_constants, in the 3 dimensions.
    """

    s = read_settings_file()['supercell_size']
    lc = a.cell.cellpar()
    return [lc[0]/s, lc[1]/s, lc[2]/s]

def volume_pressure(a):
    """Calculates volume and pressure of a material.
    Parameters:
    a (obj): a is an atoms object of class defined in ase.
    Returns:
    tuple:returns a tuple of volume and pressure.
    """
    N = len(a.get_chemical_symbols())
    vol = a.get_volume()/N
    stress = a.get_stress()
    pressure = (stress[0] + stress[1] + stress[2])/3 * units._e * units.m**3 * 10**(-9)  # eV/Å^3 to GPa
    return vol, pressure

def debye_lindemann(a, msd, temp):
    """Calculates the debye temperature and the Lindemann
       criterion. Original cell is assumed to be sc, bcc or fcc.
       Lattice constants in a,b and c may be different.
    Parameters:
    a (obj): a is an atoms object of class defined in ase.
    msd (float): mean square displacment
    temp (float): temperature
    Returns
    list: list of debye temperature and lindemann criterion.
    """
    s = read_settings_file()
    if msd != 0:
        debye = math.sqrt(9 * units._hbar**2 * temp / (units._k * a.get_masses()[0] * units._amu * msd) * units.m**2)
    else:
        debye = 0
    z = s['supercell_size']
    n = len(a) / z**3
    lc = a.cell.cellpar()[0:3]
    if n == 1:
        nnd = min(lc)
    elif n == 2:
        nnd = 1/2 * math.sqrt(lc[0]**2 + lc[1]**2 + lc[2]**2)
    elif n == 4:
        if np.max(lc) != np.min(lc): # all values in lc are not the same
            lc = np.delete(lc, np.argwhere(lc==max(lc)))
        nnd = 1/2 * math.sqrt(lc[0]**2 + lc[1]**2)
    else:
        nnd = 9999
    lindemann = math.sqrt(msd)/nnd
    return debye, lindemann

def self_diff(a, msd, time):
    """Calculates the self diffusion coefficient of a material.
    Paramters:
    a (obj): a is an atoms object of class defined in ase.
    msd (float): mean squre displacement.
    time (float): time step.
    Returns:
    float: self diffusion coefficient.
    """
    if time == 0:
        sd = "_"
    else:
        sd = msd/(6*time)
    return sd * 10 # units: mm^2 / s

def initialize_properties_file(a, ai, id, d, ma, dir=""):
    """Initializes a file over properties with correct titles and main structure
        for an material.
    Parameters:
    a (obj): a is an atoms object of class defined in ase. The material is made
            into an atoms object.
    ai (obj): initial atoms object an object of class sdefined in ase. The unit cell
                atoms object that md runs for.
    id (str): a special number identifying the material system.
    d (int): a number for the formatting of file. Give a correct spacing
            for printing to file.
    ma (boolean): a boolean indicating if the material is monoatomic
    Returns:
    None
    """
    # Help function for formating
    def lj(str, k = d):
        return " "+str.ljust(k + 6)
    file = open("property_calculations/" + dir + "properties_" + id + ".txt", "w+")

    file.write("Material ID: " + id + "\n")
    file.write("Unit cell composition: " + a.get_chemical_formula() + "\n")
    chem_formula = a.get_chemical_formula(mode='hill', empirical=True)
    file.write("Material:  "+ chem_formula + "\n")

    # Write the elements as title
    file.write("Site positions of initial unit cell:" + "\n")
    dict = chemparse.parse_formula(ai.get_chemical_formula())
    els = list(dict.keys())
    prop_num = list(dict.values())
    tmp_ls = [(a + " ") * int(b) for a,b in zip(els, prop_num)] # Get ["Al", "Mg Mg Mg"] for "AlMg3" e.g.
    els_str = "".join(tmp_ls)
    els_ls = els_str.split()  # give you ["Al", "Mg", "Mg", "Mg"] e.g.
    for a in els_ls:
        file.write(lj(a))

    # Write the site positions
    res_array = ai.get_positions()
    for i in range(0, 3): # 3 components
        file.write("\n")
        for ii in range(0, len(res_array)):
            format_str = "." + str(d) + "f"
            val  = format(res_array[:,i][ii], format_str) # d decimals
            file.write(lj(val))

    file.write("\n")
    file.write("Properties:\n")
    file.write(lj("Time")+lj("Epot")+lj("Ekin")+lj("Etot")+lj("Temp",2)+lj("MSD"))
    file.write(lj("Self_diff")+lj("LC_a",3)+lj("LC_b",3)+lj("LC_c",3))
    file.write(lj("Volume")+lj("Pressure"))
    if ma:
        file.write(lj("DebyeT",2)+lj("Lindemann")[:-1])
    file.write("\n")
    file.write(lj("fs")+lj("eV/atom")+lj("eV/atom")+lj("eV/atom")+lj("K",2)+lj("Å^2"))
    file.write(lj("mm^2/s")+lj("Å",3)+lj("Å",3)+lj("Å",3))
    file.write(lj("Å^3/atom")+lj("GPa")[:-1])
    if ma:
        file.write(lj("K",2)+lj("1"))
    file.write("\n")
    file.close()
    return

def ss(value, decimals):
    """Help function to calc_properties."""
    if isinstance(value,str):
        tmp = value
    else:
        #value = float(value)
        tmp = str(round(value, decimals))
    return " "+tmp.ljust(decimals + 6)

def calc_properties(a_old, a, id, d, ma, dir=""):
    """Calculates prioperties and writes them in a file.
    Parameters:
    a_old (obj): a_old is an atoms object from clas defined from ase.
                it's the old atom that needs to be updated.
    a (obj): a is an atoms object from clas defined from ase.
            it's the new updated atom obj for MD molecular dyanimcs.
    id (str):
d (int):
    ma (boolean):
    Returns: None
    """
    f=open("property_calculations/"+ dir +"properties_"+id+".txt", "r")
    
    #print('atoms pos', a.get_positions(), '\n')
    #print('old_atoms pos', a_old.get_positions())
    epot, ekin, etot, temp = energies_and_temp(a)
    msd =  meansquaredisp(a, a_old)
    #print('MSD', msd)
    settings = read_settings_file()
    ln = sum(1 for line in f)
    #print(ln)
    # why william?
    time = settings['time_step']*settings['interval']*(ln-10)
    #print('time', time)
    #print('time_step', settings['time_step'])
    #print('interval', settings['interval'])
    #print(ln-6)
    selfd = self_diff(a, msd, time)
    lc = lattice_constants(a)
    vol, pr = volume_pressure(a)
    f.close()

    file=open("property_calculations/"+dir+"properties_"+id+".txt", "a+")
    file.write(ss(time, d)+ss(epot, d)+ss(ekin, d)+ss(etot, d)+ss(temp, 2)+ss(msd, d))
    file.write(ss(selfd, d)+ss(lc[0], 3)+ss(lc[1], 3)+ss(lc[2], 3))
    file.write(ss(vol, 3)+ss(pr, d))
    if ma:
        debye, linde = debye_lindemann(a,msd,temp)
        file.write(ss(debye, 2)+ss(linde, d))

    file.write("\n")
    file.close()
    return 

def finalize_properties_file(a, id, d, ma, dft=False, dir="", offset=0):
    """ Calculates and records the properties of a material.
    Parameters:
    a (obj): Atoms object form ase.
    id (str): a special number identifying the material system.
    d (int): a number for the formatting of file. Give a correct appending
            for strings.
    ma (boolean): ma is a boolean, for True the system is monoatomic.
    dft (boolean): dft is True if a properties file from DFT data is to 
             be finalized.
    Returns: None
    """
    epot = []
    ekin = []
    etot = []
    temp = []
    msd = []
    selfd = []
    pr = []
    debye = []
    linde = []

    settings = read_settings_file()
    f=open("property_calculations/"+dir+"properties_"+id+".txt", "r")
    f_lines = f.readlines()
    #dft = True
    if dft:
        steps = int(8000/100)
    else:
        steps = math.floor((settings['max_steps'] - offset) / settings['interval'])
    print('properties.py: dft', dft)
    print('properties.py: steps', steps)
    for line in f_lines[-steps:]:
        #print('line', line)
        epot.append(float(line.split()[1]))
        ekin.append(float(line.split()[2]))
        etot.append(float(line.split()[3]))
        temp.append(float(line.split()[4]))
        msd.append(float(line.split()[5]))
        selfd.append(line.split()[6])
        pr.append(float(line.split()[11]))
        if ma:
            if line.split()[12] == 'inf':
                debye.append(0.0)
            else:
                debye.append(float(line.split()[12]))
            linde.append(float(line.split()[13]))
    f.close()
    
    epot_t = sum(epot)/steps
    ekin_t = sum(ekin)/steps
    etot_t = sum(etot)/steps
    temp_t = sum(temp)/steps
    msd_t = sum(msd)/steps
    try:
        selfd_t = sum(float(i) for i in selfd[1:])/(steps-1)
    except:
        selfd_t = 0
        print("selfd, division by zero")
    pr_t = sum(pr)/steps
    debye_t = sum(debye)/steps
    linde_t = sum(linde)/steps
    Cv = specific_heat(temp, len(a.get_chemical_symbols()), a)

    file=open("property_calculations/"+dir+"properties_"+id+".txt", "a+")
    file.write("\nTime averages:\n")

    # Help function for formating
    def lj(str, k = d):
        return " "+str.ljust(k+6)

    file.write(lj(" ")+lj("Epot")+lj("Ekin")+lj("Etot")+lj("Temp",2)+lj("MSD"))
    file.write(lj("Self_diff")+lj("Pressure"))
    file.write(lj("Spec_heat"))
    if ma:
        file.write(lj("DebyeT",2)+lj("Lindemann"))

    file.write("\n")

    file.write(lj(" ")+lj("eV/atom")+lj("eV/atom")+lj("eV/atom")+lj("K",2)+lj("Å^2"))
    file.write(lj("mm^2/s")+lj("GPa"))
    file.write(lj("J/(K*Kg)"))

    if ma:
        file.write(lj("K",2)+lj("1"))
    file.write("\n")

    file.write(lj(" ")+ss(epot_t, d)+ss(ekin_t, d)+ss(etot_t, d)+ss(temp_t, 2)+ss(msd_t, d))
    file.write(ss(selfd_t, d)+ss(pr_t, d))
    file.write(ss(Cv, d))
    if ma:
        file.write(ss(debye_t, 2)+ss(linde_t, d))
    file.close()
    return

def delete_properties_file(id):
    """ Deletes a property file by its id
    Parameters:
    id (): a special number identifying the material system, as an int.
    Returns: None
    """
    os.remove("property_calculations/properties_"+str(id)+".txt")
    return

def clean_property_calculations():
    """ Idea: delete all propeties files without 'Time averages:'
    in them.
    """
    print(" -- Cleaning property_calculations directory -- ")
    counter = 0
    for filename in os.listdir("property_calculations"):
        f = open("property_calculations/"+ str(filename), "r")
        if "Time averages:" not in f.read():
            counter += 1
            os.remove("property_calculations/"+str(filename))

    print(" -- Removed " + str(counter) + " properties files -- ")

def get_averaged_properties(filename, element='Al', offset=0):
    '''
    Returns averaged properties from a properties file.

    Parameters:
    filename (str): filename of the properties file to calculate 
    averages for
    offset (int): an offset for the time averaging; the time average
    starts at the offset
    '''
    # convert offset in timesteps to offest in amount of lines
    # (assuming that properties are calculated each 100 timesteps)
    offset = int(offset / 100)

    lines = []
    with open('property_calculations/'+filename, 'r') as f:
        for i, line in enumerate(f):
            if i > 10 + offset:
                if line == "Time averages:\n":
                    #print('break time')
                    break
                lines.append(line.split())
    lines = lines[:-1]
    #print('lines', len(lines))

    # create atoms object
    if element == 'Al': 
        atom = bulk('Al', 'fcc', a=4.0479, cubic=True)
    elif element == 'Si':
        atom = bulk('Si', 'fcc', a=5.4310)
    else:
        raise ValueError('Unsupported element chosen: ' + element)

    atoms = atom * 2 * (1,1,1)
    print(element, len(atoms))

    MSDs = []
    temps = []
    E_tots = []
    Cvs = []
    for line in lines:
        MSDs.append(float(line[5]))
        E_tots.append(float(line[3]))
        temps.append(float(line[4]))
        Cvs.append(specific_heat_NVT(E_tots, len(atoms), atoms, temps[-1]))
        
    MSD_averaged = []
    MSD = 0
    for i, m in enumerate(MSDs):
        if i != 0:
            MSD += m
            MSD_averaged.append(MSD / i)
        else:
            MSD_averaged.append(0)
    #print('MSD averaged', MSD_averaged)

    Cvs_averaged = []
    Cv = 0
    for i, c in enumerate(Cvs):
        Cv += c
        Cvs_averaged.append(Cv / (i+1))

    E_tots_averaged = []
    E_tot = 0
    for i, e in enumerate(E_tots):
        E_tot += e
        E_tots_averaged.append(E_tot / (i+1))
        
    return MSD_averaged, Cvs_averaged, E_tots_averaged

if __name__ == "__main__":
    #clean_property_calculations()
    M, C, E = get_averaged_properties('properties_Al_DFT_eq_2000.txt', 'Al')
    
    print(C)
    #vec_0 = np.zeros(3)
    #vec_1 = np.array([8.0803606, 8.0780096, 8.0762162])

    #v = vec_0 - vec_1
    #print(normalize_half(v))
