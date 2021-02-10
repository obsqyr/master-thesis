#!/usr/bin/env python3
from tutorial_data import compounds
from qml.kernels import gaussian_kernel
from tutorial_data import energy_pbe0
from qml.math import cho_solve
import numpy as np
from ase.io.vasp import *

#import maml
from maml.apps.pes import MTPotential

def vasp_read(directory):
    '''
    Reads from a VASP directory and produces a data set of
    positions and forces (eventually more).

    Parameters:
    directory (str): name of the targeted directory
    '''

    print("reading VASP")
    #print(read_vasp(directory + "POSCAR"))
    #print(directory + "OUTCAR")
    read_vasp_out(directory + "OUTCAR", index=0)
    #print(read_vasp_xdatcar(directory + "XDATCAR")[0])

def train_qml_regressor():
    print("training")
    
    for mol in compounds:
        mol.generate_coulomb_matrix(size=23, sorting="row-norm")
    # sin matris representation för periodiska system
        #mol.generate_fchl_representation(size=23, cut_off=10.0)
    
    # dela Al data set i train och test, beräkna MAE

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
    
if __name__ == "__main__":
    #alpha, sigma = train_qml_regressor()
    #evaluate_qml_regressor(alpha, sigma)
    #train_MTP()
    vasp_read("Al_300K/")
