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

def get_MTP_MAEs():
    dirs = [x[0] for x in os.walk('test_results')]
    print(dirs)
    dirs.remove('test_results')
    print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)

    #print(filenames)
    # this could be automated
    #filenames_al = [filenames[0][0:9], filenames[1][0:9]]
    #filenames_si = [filenames[0][9:18], filenames[1][9:18]]
    filenames_al = []
    filenames_si = []
    for filename in filenames:
        filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si.append([f for f in filename if 'Si' in f])
    print(filenames_al)
    print(filenames_si)

    # Aluminum
    data_al = []
    for filename in filenames_al:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_al.append(d)

    MAEs_al = []
    for d in data_al:
        MAEs_al.append([extract_MAE(x) for x in d])
    print(len(MAEs_al[0]))
    
    # Silicon
    data_si = []
    for filename in filenames_si:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si.append(d)

    MAEs_si = []
    for d in data_si:
        MAEs_si.append([extract_MAE(x) for x in d])
    print(len(MAEs_si[0]))

    return MAEs_al, MAEs_si


def plot_mtp_closer_to_zero():
    dirs = [x[0] for x in os.walk('test_results')]
    dirs.remove('test_results')
    #print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)

    #print(filenames)
    # this could be automated
    #filenames_al = [filenames[0][0:9], filenames[1][0:9]]
    #filenames_si = [filenames[0][9:18], filenames[1][9:18]]
    filenames_al = []
    filenames_si = []
    for filename in filenames:
        filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si.append([f for f in filename if 'Si' in f])
    #filenames_al_06 = filenames_al[1]
    #filenames_al_06 = [i for i in filenames_al_06 if not i[19:23].isdecimal()]

    filenames_al_small = []
    for filename in filenames_al:
        filenames_al_small.append([i for i in filename if not i[19:23].isdecimal()]) 
    
    #print(filenames_al)
    print(len(filenames_al_small))
    #print(filenames_si)

    filenames_si_small = []
    for filename in filenames_si:
        filenames_si_small.append([i for i in filename if not i[19:23].isdecimal()])

    #print('filenames_si', filenames_si_small)
    #filenames_al = [i[9:18] for i in filenames_al]
    #filenames_si = [i[9:18] for i in filenames_si]
    #print(filenames_al)
    #print(len(filenames_si[1]))

    # Aluminum
    data_al_small = []
    for filename in filenames_al_small:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_al_small.append(d)

    MAEs_al = []
    for d in data_al_small:
        MAEs_al.append([extract_MAE(x) for x in d])
    #print(len(MAEs_al[2]))
    #MAEs_al.append(MAEs_al[0])
    del MAEs_al[0]
    #del MAEs_al[1]
    
    # plot data_al
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(100,1000,100)

    # energy
    for pot in MAEs_al:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Al_energy_per_atom_MAE_mtp_small.png')
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
    for pot in MAEs_al:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Al_forces_MAE_mtp_small.png')
    #plt.show()

    # Silicon
    data_si_small = []
    for filename in filenames_si_small:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si_small.append(d)

    MAEs_si = []
    for d in data_si_small:
        MAEs_si.append([extract_MAE(x) for x in d])
    print('MAEs_si', MAEs_si)
    del MAEs_si[0]
    #del MAEs_al[1]
    
    # plot data_si_small
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(100,1000,100)

    # energy
    for pot in MAEs_si:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp','10.mtp'])
    plt.savefig('figures/Si_energy_per_atom_MAE_mtp_small.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    for pot in MAEs_si:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_forces_MAE_mtp_small.png')
    #plt.show()

if __name__ == "__main__":
    plot_mtp_closer_to_zero()
    
    '''
    dirs = [x[0] for x in os.walk('test_results')]
    dirs.remove('test_results')
    print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)

    #print(filenames)
    # this could be automated
    #filenames_al = [filenames[0][0:9], filenames[1][0:9]]
    #filenames_si = [filenames[0][9:18], filenames[1][9:18]]
    filenames_al = []
    filenames_si = []
    for filename in filenames:
        filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si.append([f for f in filename if 'Si' in f])
    print(filenames_al)
    print(filenames_si)

    # Aluminum
    data_al = []
    for filename in filenames_al:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_al.append(d)

    MAEs_al = []
    for d in data_al:
        MAEs_al.append([extract_MAE(x) for x in d])
    print(len(MAEs_al[0]))
    MAEs_al.append(MAEs_al[0])
    del MAEs_al[0]

    # Silicon
    data_si = []
    for filename in filenames_si:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si.append(d)

    MAEs_si = []
    for d in data_si:
        MAEs_si.append([extract_MAE(x) for x in d])
    print(len(MAEs_si[0]))

    # plot data_al
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(1000,10000,1000)

    # energy
    for pot in MAEs_al:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp', '14.mtp'])
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
    for pot in MAEs_al:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp','10.mtp', '14.mtp'])
    plt.savefig('figures/Al_forces_MAE_mtp.png')
    #plt.show()

    # plot data_si
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(1000,10000,1000)

    # no 14.mtp
    del MAEs_si[0]

    # energy
    for pot in MAEs_si:
        MAE = [m[1] for m in pot]
        #print(MAE)
        #print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)
    
    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_energy_per_atom_MAE_mtp.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    for pot in MAEs_si:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_forces_MAE_mtp.png')
    #plt.show()
    '''
