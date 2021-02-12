#!/usr/bin/env python3
from tutorial_data import compounds
from qml.kernels import gaussian_kernel
from tutorial_data import energy_pbe0
from qml.math import cho_solve
import numpy as np
from qml.representations import *

import vasp_parser as vp

#import maml
from maml.apps.pes import MTPotential

def train(X_train, Y_train):
    '''
    Trains a KRR machine for data X_train (usually positions
    representations) and target Y_train (usually potentials).

    Parameters:
    X_train (np.array): np.array of training data
    Y_train (np.array): np.array of target data

    Returns:
    alpha (np.array): returns resulting alpha values as 1D np.array
    of length len(Y_train)
    sigma (float): returns sigma used in kernel
    '''
    print("training")

    sigma = 4000
    K = gaussian_kernel(X_train, X_train, sigma)
    
    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8

    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, Y_train) 
    
    #print("Alphas:")
    #print(alpha)

    return alpha, sigma

def evaluate(alpha, sigma, X_train, X_test, Y_test):
    print("evaluating")
    Ks = gaussian_kernel(X_test, X_train, sigma)

    Y_pred = np.dot(Ks, alpha)

    print(len(Y_pred), len(Y_test))

    print("MAE:", np.mean(np.abs(Y_pred - Y_test)))

    
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
    #print(K)
    
    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8

    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, Y_training) 
    print(len(K), len(Y_training))
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
    '''
    Takes position or force data read from VASP and returns 
    it as a np.array with shape (num_atoms, timesteps, 3), where
    the 3 is for each component of the position/force.

    Parameters:
    data (np.array): np.array of shape (timesteps*num_atoms, 3)
    of positions or forces.
    num_atoms (int): number of atoms in the system.

    Returns:
    np.array: a np.array is returned of shape (num_atoms, timesteps, 3)
    '''
    divided_data = [[] for l in range(num_atoms)]

    for i in range(0, len(data), num_atoms):
        for j in range(num_atoms):
            divided_data[j].append(data[i+j])

    return np.array(divided_data)

def generate_representations(data, num_atoms, timesteps, nuclear_charge):
    '''
    Generates representation for the system at each timestep
    (currently only coulomb matrix is implemented).

    Parameters:
    data (np.array): np.array of shape (num_atoms, timesteps, 3),
    usually positions.
    num_atoms (int): number of atoms in the system
    timesteps (int): total amount of timesteps
    nuclear_charge (int): the nuclear_charge of the atoms in the system

    Returns:
    np.array: an np.array of representations

    '''
    print("generating representations")

    # generate sublists of the positions for all atoms in the 
    # system at each timestep
    positions = []
    for i in range(timesteps):
        positions.append(data[:,i])
    
    # create represenation for each timestep
    representations = []
    nuclear_charges = np.array([nuclear_charge]*num_atoms)
    for pos in positions:
        representations.append(generate_coulomb_matrix(nuclear_charges, pos, size=32, sorting="row-norm"))

    return np.array(representations)

if __name__ == "__main__":
    # import data from infiles
    forces, positions, potentials, num_atoms = vp.read_infiles('Al_300K/')

    # amount of timesteps is equal to length of potentials
    timesteps = len(potentials)
    #print(timesteps)

    # divide forces and positions into matrices for each atom
    forces = divide_data(forces, num_atoms)
    positions = divide_data(positions, num_atoms)

    # write function that converts data into representaionts
    # generate_representation(forces, 'cm')
    X_pos = generate_representations(positions, num_atoms, timesteps, 13)

    # divide into training and test data
    X_pos_train = X_pos[:250]
    X_pos_test = X_pos[250:]
    
    train_pot = potentials[:250]
    test_pot = potentials[250:]

    #print(train_f)
    #train_qml_regressor()
    alpha, sigma = train(X_pos_train, train_pot)
    evaluate(alpha, sigma, X_pos_train, X_pos_test, test_pot)
    
    # 10 ts
    # MAE: 0.053403660380666906
    # 100 timesteps
    # MAE: 2.9544048360212387
    # 200 timesteps
    # MAE: 1.2737607237465731
    # 300 ts
    # MAE: 1.3170717980190503
    # 400 ts
    # MAE: 1.3331982309170984
    # 500 ts
    # 1.3232172627262855
