#!/usr/bin/env python3

from ase.calculators.calculator import Calculator
from ase.lattice.cubic import FaceCenteredCubic
from ase.build import bulk
import numpy as np
import qml
from qml.representations import *

import train_qml as tr

class zero_calculator(Calculator):
    def get_potential_energy(self, atoms=None, force_consistent=False):
        return 0.0

    def get_forces(self, atoms):
        return np.zeros((len(atoms), 3))

    def get_stress(self, atoms):
        return np.zeros(6)

    def calculation_required(self, atoms, quantities):
        return False

class KRR_calculator(Calculator):
    def __init__(self, element, timesteps):
        print("Initializing KRR calculator")
        print("Element: " + element + ". Timesteps: " + str(timesteps))
        alphas_pot = np.loadtxt('machines/KRR/potential/' + element + '/' + str(timesteps) + '.txt')
        #print(alphas_pot)

    def get_potential_energy(self, atoms=None, force_consistent=False):
        print('atoms', atoms)
        print('atoms len', len(atoms.get_chemical_symbols()))
        N = len(atoms.get_chemical_symbols())
        X = tr.generate_sine_representation(atoms, N)
        
        return 0.0

    def get_forces(self, atoms):
        return np.zeros((len(atoms), 3))

    def get_stress(self, atoms):
        return np.zeros(6)

    def calculation_required(self, atoms, quantities):
        return False

if __name__ == "__main__":
    # Set up a crystal
    '''
    size = 5
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol="Ar",
                              size=(size, size, size),
                              pbc=True)
    '''
    atom = bulk('Al', 'fcc', a=4.0479, cubic=True)
    atoms = atom*(2,2,2)

    #print(atoms)
    x = KRR_calculator('Al', 1000)
    print(x.get_potential_energy(atoms=atoms))
