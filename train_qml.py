#!/usr/bin/env python3
from tutorial_data import compounds
from qml.kernels import gaussian_kernel
from tutorial_data import energy_pbe0
from qml.math import cho_solve
import numpy as np

def train_qml_regressor():
    print("training")
    
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

if __name__ == "__main__":
    alpha, sigma = train_qml_regressor()
    evaluate_qml_regressor(alpha, sigma)
