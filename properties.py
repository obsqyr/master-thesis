#!/usr/bin/env python3

import math
from ase import Atoms
from ase import units
import numpy as np
import os

# This file contains functions to calculate material properties

def specific_heat(temp_store, N, atoms, size):
    """Calculates the specific heat for a material.
    Given by the formula: (E[T²?] - E[T]²?)/ E[T]²? = 3*2^-1*N^-1*(1-3*N*Kb*2^-1*Cv^-1).
    Where Kb is boltzmansconstant, N is the total number of atoms, T is temperature and Cv the specific heat.
    E[A(t)] calculates the expectation value of A, which can in this case be seen as a time average for the
    phase variable A(t).
    Parameters:
    temp_store (list): The list over all intantaneous temperatures of a material once MD has been run.
    N (int): The total number of atoms in the material.
    Returns:
    float: specific heat is returned (J/(K*Kg))
    """
    
    if len(temp_store) == 0:
        raise ValueError("temp_store is empty, invalid value.")
    steps = len(temp_store)
    z = sum(atoms.get_masses()) * units._amu # total mass: atomic units to kg
    ET = sum(temp_store)/steps
    ET2 = sum(np.array(temp_store)**2)/steps
    N = N / size**3
    Cv = ((9*ET**2*N*units._k) / (ET**2 * (6+4*N) - 4*N*ET2)) / z*size**3
    return Cv

