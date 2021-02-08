#!/usr/bin/env python3

from ase.calculators.calculator import Calculator
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
import qml
from qml.representations import *

alpha = np.loadtxt("machine.txt")

class qml_calculator(Calculator):
    
    def get_potential_energy(self, atoms=None, force_consistent=False):
        coordinates = atoms.get_positions()
        #print(len(coordinates))
        nuclear_charges = np.ones(len(coordinates)) 
        #print(atoms.get_charges())
        #print(nuclear_charges)

        #cm = generate_coulomb_matrix(nuclear_charges, coordinates,
        #                             size=5, sorting="row-norm")

        #print(cm)
        #Ks = gaussian_kernel(X_test, X_training, sigma)

        # Make the predictions
        #Y_predicted = np.dot(Ks, alpha)

        return 0.0

    def get_forces(self, atoms):
        return np.zeros((len(atoms), 3))

    def get_stress(self, atoms):
        return np.zeros(6)

    def calculation_required(self, atoms, quantities):
        return False

if __name__ == "__main__":
    # Set up a crystal
    size = 5
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol="Ar",
                              size=(size, size, size),
                              pbc=True)
    #print(atoms)
    x = qml_calculator()
    print(x.get_potential_energy(atoms=atoms))
