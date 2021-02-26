#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import os

def extract_MAE(data):
    '''
    Extracts all Average absolute differences / MAEs
    from a test results file.
    
    Input:
    data (str): test results file as string

    Returns:
    list: returns list of floats of all 5 MAEs.
    Index 0: Energy
    Index 1: Energy per atom
    Index 2: Forces
    Index 3: Stresses (in eV)
    Index 4: Stresses (in GPa)
    '''
    MAE = []
    for line in data.split('\n'):
        if "Average absolute difference" in line:
            MAE.append(float(line.strip().split('=')[1][1:]))
    return MAE


if __name__ == "__main__":
    dirs = [x[0] for x in os.walk('test_results')]
    dirs.remove('test_results')
    print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)

    print(filenames)

    data = []
    for filename in filenames:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data.append(d)

    MAEs = []
    for d in data:
        MAEs.append([extract_MAE(x) for x in d])
    print(len(MAEs[0]))

    
    # plot data
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(1000,10000,1000)

    # energy
    for pot in MAEs:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Al_energy_per_atom_MAE_mtp.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    for pot in MAEs:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp','10.mtp'])
    plt.savefig('figures/Al_forces_MAE_mtp.png')
    #plt.show()
    
    '''dir = "test_results/06/"
    filenames = [f for f in os.listdir(dir)]
    filenames = sorted(filenames)
    
    data = []
    for filename in filenames:
        with open(dir+filename, 'r') as f:
            #print(f.read())
            data.append(f.read())

    MAEs = [extract_MAE(d) for d in data]
    
    #print(len(MAEs))

    # plot data
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al, 06.mtp")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(1000,10000,1000)

    # energy
    MAE = [m[1] for m in MAEs]
    #print(MAE)
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    plt.savefig('figures/Al_energy_per_atom_MAE_mtp.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, Al, 06.mtp")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    
    MAE = []
    MAE = [m[2] for m in MAEs]
    #print(MAE)
    
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    plt.savefig('figures/Al_forces_MAE_mtp.png')
    #plt.show()

    #print(extract_MAE(data))

    #print(data)
'''
