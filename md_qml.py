#!/usr/bin/env python3
"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase import units
from ase.build import bulk
from asap3 import LennardJones
from asap3 import EMT
from asap3 import Trajectory
import qml
import numpy as np
import copy
import time

import properties as pr
import calculators as calcs
from read_settings import read_settings_file

def generate_qml_potential():
    alpha = np.loadtxt("machine.txt")

def run_md(calculator, timesteps, element='Al', mtp='06', eq=0, dir=""):
    # Read settings
    settings = read_settings_file('settings.json')
    # Scale atoms object, cubic
    size = settings['supercell_size']
    if element == 'Al':
        atom = bulk('Al', 'fcc', a=4.0479, cubic=True)
    elif element == 'Si':
        atom = bulk('Si', 'fcc', a=5.4310)
    else:
        raise ValueError("Invalid element chosen")
        
    initial_unitcell_atoms = copy.deepcopy(atom)
    atoms = atom * size * (1,1,1)

    # Set up a crystal
    # Al crystal, 32 atoms
    #atoms = atom*(2,2,2)
    #size = 10
    #atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Al", size=(size, size, size), pbc=True)
    N = len(atoms.get_chemical_symbols())
    print('N', N)
    
    # Create a copy of the initial atoms object for future reference
    old_atoms = copy.deepcopy(atoms)

    # Describe the interatomic interactions with QML
    #atoms.calc = calcs.zero_calculator()

    if calculator == 'MTP':
        atoms.calc = calcs.MTP_calculator(element, timesteps, mtp, eq)
    elif calculator == 'KRR':
        atoms.calc = calcs.KRR_calculator(element, timesteps)
    else:
        raise ValueError('Invalid calculator chosen')

    #atoms.calc = calcs.zero_calculator()
    # Set the momenta corresponding to T=300K
    print('atoms', atoms)
    # do I need this?
    MaxwellBoltzmannDistribution(atoms, temperature_K = 300)

    # Select integrator
    if settings['ensemble'] == "NVE":
        from ase.md.verlet import VelocityVerlet
        dyn = VelocityVerlet(atoms, settings['time_step'] * units.fs)
        
    elif settings['ensemble'] == "NVT":
        print('-- Using Langevin, NVT ensamble --')
        from ase.md.langevin import Langevin
        dyn = Langevin(atoms, settings['time_step'] * units.fs, settings['temperature'] * units.kB, settings['friction'])
        #from ase.md.nptberendsen import NPTBerendsen
        #dyn = NVTBerendsen(atoms, 1 * units.fs, 300, taut=0.5*1000*units.fs)
        #dyn = NPTBerendsen(atoms, timestep=0.1 * units.fs, temperature_K=300, taut=100 * units.fs, pressure_au=1.01325 * units.bar, taup=1000 * units.fs, compressibility=4.57e-5 / units.bar)
        
        #from ase.md.andersen import Andersen
        #dyn = Andersen(atoms, settings['time_step'] * units.fs, settings['temperature'], 0.1)
        #print("-- Using Nose-Hoover, 'NTP' ensamble --")
        #dyn = NPT(atoms, settings['time_step'] * units.fs, settings['temperature'] * units.kB, externalstress = 0, ttime = 25*units.fs, pfactor=(400*units.fs)**2 * 0.6) 
        
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    #dyn = VelocityVerlet(atoms, settings['time_step'] * units.fs)  # 5 fs time step.
    #dyn.attach(traj.write, interval=10)

    # Number of decimals for most calculated properties.
    decimals = settings['decimals']
    # Boolean indicating if the material is monoatomic.
    monoatomic = len(set(atoms.get_chemical_symbols())) == 1

    # Calculation and writing of properties
    # element + calculator as id
    if calculator == 'MTP':
        id = element + '_' + calculator + '_' + mtp + '_' + str(timesteps) +'_eq_' + str(eq) + '_ranfor_' + str(settings['max_steps'])
    elif calculator == 'KRR':
        id = element + '_' + calculator + '_' + str(timesteps)
    pr.initialize_properties_file(atoms, initial_unitcell_atoms, id, decimals, monoatomic, dir)
    dyn.attach(pr.calc_properties, 100, old_atoms, atoms, id, decimals, monoatomic, dir)

    # attach trajectory
    #traj = Trajectory('trajectories/' + id + '.traj', 'w', atoms)
    #dyn.attach(traj.write, interval=100)
    
    temperatures = []
    def printenergy(t=temperatures, a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        #print(atoms)
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        t.append(ekin / (1.5 * units.kB))
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    
    # Now run the dynamics
    dyn.attach(printenergy, interval=10)
    #dyn.attach(printenergy, interval=settings['interval'])
    #if printenergy() != None:
    #print(printenergy())

    dyn.run(settings['max_steps'])
    pr.finalize_properties_file(atoms, id, decimals, monoatomic, pr, dir)
    print('id', id)
    #return temperatures, N, atoms, size
    return id

if __name__ == "__main__":
    start_time = time.time()
    run_md('MTP', 100, 'Si', '06', 0)
    print("Ran in %s seconds" % (time.time() - start_time))
    ## calculate specific heat
    #spec_heat = properties.specific_heat(temperatures, N, atoms, size) / 1000 # convert to KJ/K*kg
    #print ("Specific heat " + str(atoms.symbols) + ": %.4f [kJ/(K*kg)]" % (spec_heat))
