#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import os
import train_qml as tr

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
    dirs = [x[0] for x in os.walk('MTP/test_results')]
    dirs.remove('MTP/test_results')
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
    #print(filenames_al)
    #print(filenames_si)

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

def scatter_plot(x, y, filename, title='', xlabel='', ylabel='', legend=''):
    print("scatter plot")
    
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    for a,b in zip(x,y):
        plt.scatter(a, b)
        plt.plot(a, b)

    plt.legend(legend)
    plt.savefig(filename)
    plt.show()

def plot_energies():
    print("plotting energies")
    MAEs_al, MAEs_si = get_MTP_MAEs()
    
    # remove incomplete 14.mtp from both data sets
    del MAEs_si[0]
    del MAEs_al[0]

    # plot Al energy
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
        
    qml_MAEs_al = np.loadtxt('potentials_MAEs_al.txt')
    plt.scatter(timesteps, qml_MAEs_al)
    plt.plot(timesteps, qml_MAEs_al)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])
    plt.savefig('figures/Al_energy_per_atom_MAE_all.png')
    #plt.show()

    # plot Si energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(1000,10000,1000)

    # energy
    for pot in MAEs_si:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    qml_MAEs_si = np.loadtxt('potentials_MAEs_Si.txt')
    plt.scatter(timesteps, qml_MAEs_si)
    plt.plot(timesteps, qml_MAEs_si)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])

    plt.savefig('figures/Si_energy_per_atom_MAE_all.png')
    #plt.show()

    
def plot_forces():
    #forces_MAEs = np.loadtxt("forces_MAEs.txt")
    #print(forces_MAEs)
    #for i in range(32):
    #    plt.scatter(i, forces_MAEs[i][0])
    
    #plt.savefig('figures/forces_x_component.png')
    
    filenames = [f for f in os.listdir("forces_MAEs")]
    filenames = sorted(filenames)
    MAEs = []
    for f in filenames:
        MAEs.append(np.loadtxt("forces_MAEs/"+f))
        
    MAEs = np.array(MAEs)
    print(MAEs.shape)

    plt.title("All atoms, Al, force lengths MAE against timesteps")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    timesteps = range(1000,10000,1000)
    
    for i in range(32):
    #print(MAEs[:,0])
        plt.scatter(timesteps, MAEs[:,i])
        plt.plot(timesteps, MAEs[:,i])
    
    #plt.legend(['x', 'y', 'z'])
    
    plt.savefig('figures/Al_force_length_MAE.png')
    plt.show()

if __name__ == "__main__":
    #plot_forces()
    plot_energies()


    # håll separata tränings- och evalueringsdataset
    # säg håll ett test

    # träna explicit för simulering, mata in säg 10 observerade värden?
    # när vi tränar, optimera för att mata in modellens egna prediktioner
    # och komma så nära det riktiga datat som möjligt
    # formulera om som ett tidsserie problem? autoregressiv modell?
    # börja med enkel KRR, jämför med detta? parallelt jämföra med MTP
    # också?
