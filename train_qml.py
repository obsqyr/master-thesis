#!/usr/bin/env python3
from tutorial_data import compounds
from qml.kernels import gaussian_kernel
from tutorial_data import energy_pbe0
from qml.math import cho_solve
import numpy as np

import vasp_parser as vp

#import maml
from maml.apps.pes import MTPotential

def train(tr, target):
    print("training")

    sigma = 4000
    K = gaussian_kernel(tr, tr, sigma)

    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8

    print(len(K), len(target))

    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, target) 

    print("Alphas:")
    print(alpha)
    
def train_qml_regressor():
    print("training QML regressor")
    
    for mol in compounds:
        mol.generate_coulomb_matrix(size=23, sorting="row-norm")
    # sin matris representation för periodiska system
        #mol.generate_fchl_representation(size=23, cut_off=10.0)
    
    # dela Al data set i train och test, beräkna MAE

    # Make a big 2D array with all the representations
    X = np.array([mol.representation for mol in compounds])
    print(len(X[0]))
    #print(X)

    # Assign 1000 first molecules to the training set
    X_training = X[:1000]
    Y_training = energy_pbe0[:1000]
   
    sigma = 4000
    K = gaussian_kernel(X_training, X_training, sigma)
    #print("Gaussian kernel:")
    
    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8

    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, Y_training) 
    #print(len(K), len(Y_training))
    print(Y_training)

    #print("Alphas:")
    #print(alpha)
    # write resulting alphas to file
    np.savetxt("machine.txt", X = alpha, header = "representation: CM, regressor: KRR")
    return alpha, sigma

def evaluate_qml_regressor(alpha, sigma):
    print("evaluating")
    
    for mol in compounds:
        mol.generate_coulomb_matrix(size=23, sorting="row-norm")
    
    # Make a big 2D array with all the representations
    X = np.array([mol.representation for mol in compounds])

    # Assign 1000 last molecules to the test set
    X_test = X[-1000:]
    Y_test = energy_pbe0[-1000:]

    # Assign 1000 first molecules to the training set
    X_training = X[:1000]
    Y_training = energy_pbe0[:1000]

    # calculate a kernel matrix between test and training data, using the same sigma
    Ks = gaussian_kernel(X_test, X_training, sigma)
    print(X_test[0])
    # Make the predictions
    Y_predicted = np.dot(Ks, alpha)

    # Calculate mean-absolute-error (MAE):
    print (np.mean(np.abs(Y_predicted - Y_test)))

def train_MTP():
    print("Training MTP")
    mtp = MTPotential()
    # kör train, använd evaluate i ett givet tidsteg för att få
    # ut värden på krafter och energi, implementera en ASE calculator
    # med detta.

def divide_data(data, num_atoms):
    divided_data = [[] for l in range(num_atoms)]

    for i in range(0, len(data), num_atoms):
        for j in range(num_atoms):
            divided_data[j].append(data[i+j])

    return np.array(divided_data)

if __name__ == "__main__":
    # import data from infiles
    forces, positions, potentials = vp.read_infiles()

    # amount of timesteps is equal to length of potentials
    timesteps = len(potentials)
    print(timesteps)

    # divide forces and positions into matrices for each atom
    forces = divide_data(forces, 32)
    positions = divide_data(positions, 32)
    
    # divide into training and test data
    train_f = forces[:,0:50]
    test_f = forces[:, 50:]
    
    train_pos = positions[:,0:50]
    test_pos = positions[:, 50:]
    
    train_pot = potentials[:50]
    test_pot = potentials[50:]

    #print(train_pot)
    #print(test_pot)
    
    #train_qml_regressor()
    #train(train_pos, train_f)
