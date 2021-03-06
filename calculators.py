#!/usr/bin/env python3

from ase.calculators.calculator import Calculator
from ase.lattice.cubic import FaceCenteredCubic
from ase.build import bulk
import numpy as np
import qml
from qml.representations import *
import subprocess
from os import path

# -- own modules --
import train_qml as tr
import MTP.cfg_parser as cfg_parser

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
        self.alphas_pot = np.loadtxt('machines/KRR/potential/' + element + '/alpha/' + str(timesteps) + '.txt')
        self.X_train = np.loadtxt('machines/KRR/potential/' + element + '/training_data/' + str(timesteps) + '.txt')
        self.alphas_forces = np.load('machines/KRR/forces/' + element + '/alpha/' + str(timesteps) + '.npy')

    def get_potential_energy(self, atoms=None, force_consistent=False):
        #print('atoms', atoms)
        #print('atoms len', len(atoms.get_chemical_symbols()))
        N = len(atoms.get_chemical_symbols())
        X = tr.generate_sine_representation(atoms, N)
        X = np.array([X])
        
        pred = tr.predict_potential(self.alphas_pot, 4000, self.X_train, X)
        return float(pred)

    def get_forces(self, atoms):
        N = len(atoms.get_chemical_symbols())
        X = tr.generate_sine_representation(atoms, N)
        X = np.array([X])
        pred = tr.predict_forces(self.alphas_forces, 4000, self.X_train, X)
        #print(pred.shape)
        return pred
        #return np.zeros((len(atoms), 3))

    def get_stress(self, atoms):
        return np.zeros(6)

    def calculation_required(self, atoms, quantities):
        return False

class MTP_calculator(Calculator):
    def __init__(self, element, timesteps, mtp, eq=0):
        #print("Initializing MTP calculator")
        #super(MTP_calculator, self).__init__()
        #self.implemented_properties = ['energy', 'forces']
        if eq == 0:
            print("Element: " + element + ". Timesteps: " + str(timesteps) + ". Potential: " + mtp)
            self.mtp_path = 'MTP/mtps_out/' + element + '_' + mtp + '_pot_' + str(timesteps) + '.mtp'
        else:
            print("Element: " + element + ". Timesteps: " + str(timesteps) + ". Potential: " + mtp + ". Eq. after: " + str(eq))
            self.mtp_path = 'MTP/mtps_out_eq/' + element + '_' + mtp + '_pot_' + str(timesteps) + '_eq_' + str(eq) + '.mtp'

        if not path.exists(self.mtp_path):
            raise ValueError('Chosen potential does not exist')
        
    def get_potential_energy(self, atoms=None, force_consistent=False):
        cfg_parser.atoms_to_cfg(atoms, 'MTP/atom.cfg')
        subprocess.call(['MTP/predict_mtp.sh', self.mtp_path])
        pot = cfg_parser.read_potential_cfg('MTP/pred.cfg')
        return pot
        
    def get_forces(self, atoms):
        cfg_parser.atoms_to_cfg(atoms, 'MTP/atom.cfg')
        subprocess.call(['MTP/predict_mtp.sh', self.mtp_path])
        N = len(atoms)
        forces = cfg_parser.read_forces_cfg('MTP/pred.cfg', N)

        return forces
        #return np.zeros((len(atoms), 3))

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
    x = MTP_calculator('Al', 1000, '06')
    print(x.get_forces(atoms=atoms))
