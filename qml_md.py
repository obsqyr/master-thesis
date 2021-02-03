#!/usr/bin/env python3
"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import LennardJones
from asap3 import EMT
from asap3 import Trajectory
import properties as pr
import qml
import numpy as np

def generate_qml_potential():
    from tutorial_data import compounds
    from qml.kernels import gaussian_kernel
    from tutorial_data import energy_pbe0
    from qml.math import cho_solve

    for mol in compounds:
        mol.generate_coulomb_matrix(size=23, sorting="row-norm")
        #mol.generate_fchl_representation(size=23, cut_off=10.0)
    
    # Make a big 2D array with all the representations
    X = np.array([mol.representation for mol in compounds])
    #print(X)

    # Assign 1000 first molecules to the training set
    X_training = X[:1000]
    Y_training = energy_pbe0[:1000]
   
    sigma = 4000
    K = gaussian_kernel(X_training, X_training, sigma)
    #print("Gaussian kernel:")
    #print(K)

    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8

    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, Y_training) 

    #print("Alphas:")
    print(alpha)

    # Assign 1000 last molecules to the test set
    X_test = X[-1000:]
    Y_test = energy_pbe0[-1000:]

    # calculate a kernel matrix between test and training data, using the same sigma
    Ks = gaussian_kernel(X_test, X_training, sigma)

    # Make the predictions
    Y_predicted = np.dot(Ks, alpha)

    # Calculate mean-absolute-error (MAE):
    print (np.mean(np.abs(Y_predicted - Y_test)))

# Set up a crystal
size = 10
atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Ar",
                          size=(size, size, size),
                          pbc=True)
N = len(atoms.get_chemical_symbols())

# Describe the interatomic interactions with the Effective Medium Theory
#atoms.calc = EMT()
# Describe the interatomic interactions with Lennard Jones potential
atoms.calc = LennardJones([18], [0.010323], [3.40], rCut = 6.625, modified = True)
# Describe the interatomic interactions with QML
generate_qml_potential()

'''
# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, temperature_K = 300)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
traj = Trajectory(atoms.get_chemical_symbols()[0] + '.traj', 'w', atoms)
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
dyn.attach(traj.write, interval=10)

temperatures = []
def printenergy(t=temperatures, a=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    t.append(ekin / (1.5 * units.kB))
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# Now run the dynamics
dyn.attach(printenergy, interval=10)
print(printenergy())
dyn.run(2000)

# calculate specific heat
spec_heat = pr.specific_heat(temperatures, N, atoms, size) / 1000 # convert to KJ/K*kg
print ("Specific heat " + str(atoms.symbols) + ": %.4f [kJ/(K*kg)]" % (spec_heat))
'''
