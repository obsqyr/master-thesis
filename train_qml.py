#!/usr/bin/env python3
#from tutorial_data import compounds
from qml.kernels import gaussian_kernel
#from tutorial_data import energy_pbe0
from qml.math import cho_solve
import numpy as np
from qml.representations import *
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory

from dscribe.descriptors import SineMatrix

from ase.build import bulk
from ase.visualize import view

# own modules
import vasp_parser as vp
import visualization as vis

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

    print(X_train.shape, Y_train.shape)
    # Use the built-in Cholesky-decomposition to solve
    alpha = cho_solve(K, Y_train) 
    
    #print("Alphas:")
    #print(alpha)

    return alpha, sigma

def evaluate(alpha, sigma, X_train, X_test, Y_test):
    print("evaluating")

    Ks = gaussian_kernel(X_test, X_train, sigma)

    Y_pred = np.dot(Ks, alpha)

    #print(len(Y_pred), len(Y_test))

    MAE = np.mean(np.abs(Y_pred - Y_test))
    print("MAE:", MAE)
    return MAE    

def train_forces(X_train, Y_train):
    print("training forces machines")

    sigma = 4000
    K = gaussian_kernel(X_train, X_train, sigma)
    
    # Add a small lambda to the diagonal of the kernel matrix
    K[np.diag_indices_from(K)] += 1e-8
    
    alphas = []
    for atom in Y_train:
        print(atom.shape)
        alphas_comp = []
        # go through all three force components
        for component in range(0,3):
            y_train = atom[:,component]
            print(X_train.shape, y_train.shape)
            # Use the built-in Cholesky-decomposition to solve
            alpha_comp = cho_solve(K, y_train) 
            alphas_comp.append(alpha_comp)
        alphas.append(alphas_comp)

    print('alphas: ', np.array(alphas).shape)
    #print("Alphas:")
    #print(alpha)

    return np.array(alphas), sigma

def evaluate_forces(alphas, sigma, X_train, X_test, Y_test):
    print("evaluating forces machines")
    
    Ks = gaussian_kernel(X_test, X_train, sigma)

    component_MAEs = []
    print(Y_test.shape)
    for i, atom in enumerate(Y_test):
        print('atom', atom.shape)
        Y_pred = []
        y_test = []
        for component in range(0,3):
            print(Ks.shape, alphas[i].shape)
            Y_pred_component = np.dot(Ks, alphas[i, component,:])
            Y_pred.append(Y_pred_component)
            
            y_test_component = atom[:,component]
            y_test.append(y_test_component)

        Y_pred = np.linalg.norm(np.array(Y_pred), axis=0)
        y_test = np.linalg.norm(np.array(y_test), axis=0)
        #print('before norm', Y_pred.shape, y_test.shape)
        #Y_pred = np.linalg.norm(Y_pred, axis=0)
        #y_test = np.linalg.norm(y_test, axis=0)
        #print('after norm', Y_pred.shape, y_test.shape)
        #print(Y_pred, y_test)

        MAE = np.mean(np.abs(Y_pred - y_test))
        print("MAE:", MAE)
        component_MAEs.append(MAE)

    return np.array(component_MAEs)
    
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
    #print(Y_training)

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

def generate_representations(data, timesteps, atoms, num_atoms, atomic_num, rep_type = 'sine'):
    '''
    Generates representation for the system at each timestep
    (currently only coulomb and sine matrices are implemented).

    Parameters:
    data (np.array): np.array of shape (num_atoms, timesteps, 3),
    usually positions.
    timesteps (int): total amount of timesteps
    num_atoms (int): number of atoms in the system
    atomic_num (int): the atomic number of the atoms in the system

    Returns:
    np.array: an np.array of representations

    '''
    print("generating " + rep_type +" matrix representations")
    
    # generate sublists of the positions for all atoms in the 
    # system at each timestep
    positions = []
    for i in range(timesteps):
        positions.append(data[:,i])
    
    if rep_type == 'cm':
        # create represenation for each timestep
        representations = []
        print(type(num_atoms))
        atomic_numbers = np.array([atomic_num]*num_atoms)
        #print(nuclear_charges)
        # högre ordnings, mata in säg 10 tidsteg, prediktera
        for pos in positions:
            representations.append(generate_coulomb_matrix(atomic_numbers, pos, size=32, sorting="row-norm"))
            #representations.append(SineMatrix(n_atoms_max=num_atoms, permutation='none', sparse=False, flatten=True))
            #print(representations[0])
        return np.array(representations)
    elif rep_type == 'sine':    
        # create sine matrix descriptor for system
        sm = SineMatrix(n_atoms_max=num_atoms, permutation='none', sparse=False, flatten=True)
        X = [sm.create(x) for x in atoms]
        return np.squeeze(np.array(X))
    
def generate_sine_matrix(data, timesteps, num_atoms):
    print("generating sine matrix representation")

    # generate sublists of the positions for all atoms in the 
    # system at each timestep
    positions = []
    for i in range(timesteps):
        positions.append(data[:,i])

    # hur vill jag göra det här? skriv ner atoms objekt i read vasp
    # eller skapa nytt atoms objekt som nedan?
    system = []
    atom = bulk('Al', 'fcc', a=4.0479, cubic=True)
    atoms = atom*(2,2,2)

    # generate list of atoms objects with positions for each timestep
    for pos in positions:
        atoms.set_positions(pos)
        system.append(atoms.copy())
    
    # create sine matrix descriptor for system
    sm = SineMatrix(n_atoms_max=num_atoms, permutation='none', sparse=False, flatten=True)
    X = [sm.create(x) for x in system]
    return np.squeeze(np.array(X))

def train_and_plot_potentials_machine(X_pos, potentials):
    # reserve last 1000 data points for testing
    test_pot = potentials[9000:]
    X_pos_test = X_pos[9000:]
    
    #print(train_f)
    #train_qml_regressor()
    indeces = range(0,10000, 1000)
    print(indeces)
    MAEs = []
    for i in indeces[1:]:
        X_pos_train = X_pos[:i]
        train_pot = potentials[:i]
        print(len(X_pos_train), len(train_pot))
        
        alpha, sigma = train(X_pos_train, train_pot)
        MAEs.append(evaluate(alpha, sigma, X_pos_train, X_pos_test, test_pot))

    print(MAEs)
    MAEs = [x / 32 for x in MAEs]
    x = []
    y = []
    x.append(indeces[1:])
    y.append(MAEs)
    
    # reserve last 8000 data points for testing
    test_pot = potentials[8000:]
    X_pos_test = X_pos[8000:]
    
    #print(train_f)
    #train_qml_regressor()
    indeces = range(0,9000,1000)
    MAEs = []
    for i in indeces[1:]:
        X_pos_train = X_pos[:i]
        train_pot = potentials[:i]
        print(len(X_pos_train), len(train_pot))
        
        alpha, sigma = train(X_pos_train, train_pot)
        MAEs.append(evaluate(alpha, sigma, X_pos_train, X_pos_test, test_pot))

    print(MAEs)
    MAEs = [x / 32 for x in MAEs]
    x.append(indeces[1:])
    y.append(MAEs)

    vis.scatter_plot(x, y, "figures/MAE_1000and2000points_Al.png", "MAE against timesteps (sine matrix, Al)", 'timesteps', 'MAE [eV/atom]', ['1000 test data points', '2000 test data points'])

def train_and_evaluate_forces(X_pos, forces, indeces):
    for i in indeces[1:]:
        alphas, sigma = train_forces(X_pos[:i], forces[:,:i])
        MAEs = evaluate_forces(alphas, sigma, X_pos[:i], X_pos[9000:], forces[:,9000:])    
        print(MAEs.shape)
        print(MAEs)
        np.savetxt("forces_MAEs/" + str(i) + ".txt", MAEs)

if __name__ == "__main__":
    # import data from infiles
    forces, positions, potentials, atom_prop = vp.read_infiles('Al_300K/')

    # convert np.arrays to ints
    num_atoms = int(atom_prop[0])
    atomic_num = int(atom_prop[1])
    
    # import atoms from trajectory file
    traj = Trajectory('Al_300K_infiles/Al.traj')
    atoms = [atom for atom in traj]
    
    # amount of timesteps is equal to length of potentials
    timesteps = len(potentials)

    # divide forces and positions into matrices for each atom
    forces = divide_data(forces, num_atoms)
    positions = divide_data(positions, num_atoms)

    X_pos = generate_representations(positions, timesteps, atoms, num_atoms, atomic_num, 'sine')

    #print(np.linalg.norm(forces[0]))

    #train(X_pos, potentials[:500])
    train_and_evaluate_forces(X_pos, forces, range(0,10000,1000))
    
    # consider randomising, evaluate on other data
    # time series, transformer model? samla in mycket data
    # divide into training and test data
    #X_pos_train = X_pos[:5000]
    #X_pos_test = X_pos[5000:]
    
    #train_pot = potentials[:250]

    #train_and_plot_potentials_machine(X_pos, potentials)
