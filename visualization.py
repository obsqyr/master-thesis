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

def get_MTP_MAEs(size='big'):
    '''
    size: 
        'big' - MAEs for timesteps 1000, 2000, ... , 9000
        'small' - MAEs for timesteps 100, 200, ... , 900
        'smaller' - MAEs for timesteps 10, 20, ... , 90
    '''
    dirs = [x[0] for x in os.walk('MTP/test_results')]
    dirs.remove('MTP/test_results')
    #print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)

    # this could be automated
    #filenames_al = [filenames[0][0:9], filenames[1][0:9]]
    #filenames_si = [filenames[0][9:18], filenames[1][9:18]]
    # extract all filenames
    filenames_al = []
    filenames_si = []
    for filename in filenames:
        filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si.append([f for f in filename if 'Si' in f])
    #print('filenames_al', filenames_al)
    #print(filenames_si)

    # resize all filenames to chosen size
    if size == 'big':
        temp_al = []
        for filename in filenames_al:
            #print(filename)
            temp_al.append([i for i in filename if i[23:27].isdecimal()])
        filenames_al = temp_al
        temp_si = []
        for filename in filenames_si:
            #print(filename)
            temp_si.append([i for i in filename if i[23:27].isdecimal()])
        filenames_si = temp_al
    elif size == 'small':
        temp_al = []
        for filename in filenames_al:
            #print(filename)
            temp_al.append([i for i in filename if i[23:26].isdecimal() and not i[23:27].isdecimal()])
        filenames_al = temp_al
        temp_si = []
        for filename in filenames_si:
            #print(filename)
            temp_si.append([i for i in filename if i[23:26].isdecimal() and not i[23:27].isdecimal()])
        filenames_si = temp_si
    elif size == 'smaller':
        temp_al = []
        for filename in filenames_al:
            #print(filename)
            temp_al.append([i for i in filename if i[23:25].isdecimal() and not i[23:26].isdecimal()])
        filenames_al = temp_al
        temp_si = []
        for filename in filenames_si:
            #print(filename)
            temp_si.append([i for i in filename if i[23:25].isdecimal() and not i[23:26].isdecimal()])
        filenames_si = temp_si

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

def plot_energies(size='big'):
    print("plotting energies")
    MAEs_al, MAEs_si = get_MTP_MAEs(size)
    
    # remove incomplete 14.mtp from both data sets
    del MAEs_si[0]
    del MAEs_al[0]

    # plot Al energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    # this depends on size
    if size == 'big':
        timesteps = range(1000,10000,1000)
    elif size == 'small':
        timesteps = range(100, 1000, 100)
    elif size == 'smaller':
        timesteps = range(10, 100, 10)

    # energy
    for pot in MAEs_al:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)
        
    qml_MAEs_al = np.loadtxt('potentials_MAEs/Al_' + size + '.txt')
    plt.scatter(timesteps, qml_MAEs_al)
    plt.plot(timesteps, qml_MAEs_al)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])
    plt.savefig('figures/Al_energy_per_atom_MAE_all_' + size + '.png')
    plt.show()

    # plot Si energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')

    # energy
    for pot in MAEs_si:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    qml_MAEs_si = np.loadtxt('potentials_MAEs/Si_' + size + '.txt')
    plt.scatter(timesteps, qml_MAEs_si)
    plt.plot(timesteps, qml_MAEs_si)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])

    plt.savefig('figures/Si_energy_per_atom_MAE_all_' + size + '.png')
    #plt.show()

def plot_forces(size='big'):
    print("plotting forces")
    MAEs_al, MAEs_si = get_MTP_MAEs(size)
    
    # remove incomplete 14.mtp from both data sets
    del MAEs_si[0]
    del MAEs_al[0]

    # plot Al forces
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Average force length, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/ Å]')
    # this depends on size
    if size == 'big':
        timesteps = range(1000,10000,1000)
    elif size == 'small':
        timesteps = range(100, 1000, 100)
    elif size == 'smaller':
        timesteps = range(10, 100, 10)

    # forces
    for pot in MAEs_al:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)
        
    qml_MAEs_al = np.loadtxt('forces_MAEs/Al_' + size + '.txt')
    plt.scatter(timesteps, qml_MAEs_al)
    plt.plot(timesteps, qml_MAEs_al)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])
    #plt.grid(True)
    plt.savefig('figures/Al_mean_forces_lengths_MAE_all_' + size + '.png')
    plt.show()

    # plot Si energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Average force length, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV / Å]')

    # forces
    for pot in MAEs_si:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    qml_MAEs_si = np.loadtxt('forces_MAEs/Si_' + size + '.txt')
    plt.scatter(timesteps, qml_MAEs_si)
    plt.plot(timesteps, qml_MAEs_si)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])

    plt.savefig('figures/Si_mean_forces_lengths_MAE_all_' + size + '.png')
    #plt.show()

def plot_forces_sphere():
    #forces_MAEs = np.loadtxt("forces_MAEs.txt")
    #print(forces_MAEs)
    #for i in range(32):
    #    plt.scatter(i, forces_MAEs[i][0])
    
    #plt.savefig('figures/forces_x_component.png')
    
    filenames = [f for f in os.listdir("forces_MAEs")]
    filenames = sorted(filenames)
    filenames_sphere = filenames[9:18]
    filenames = filenames[0:9]
    print(filenames)
    print(filenames_sphere)
    MAEs = []
    for f in filenames:
        MAEs.append(np.loadtxt("forces_MAEs/"+f))
        
    MAEs_s = []
    for f in filenames_sphere:
        MAEs_s.append(np.loadtxt("forces_MAEs/"+f))
    
    MAEs = np.array(MAEs)
    MAEs_s = np.array(MAEs_s)
    print('MAEs', MAEs.shape, 'MAEs_s', MAEs_s.shape)

    plt.title("Mean over all atoms, Al, force length MAE against timesteps")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    timesteps = range(1000,10000,1000)
    
    MAEs = np.mean(MAEs, axis=1)
    MAEs_s = np.mean(MAEs_s, axis=1)

    print('mean MAEs', MAEs.shape)
    
    # for i in range(32):
    # plt.scatter(timesteps, MAEs[:,i], color='r')
    # etc ...
    plt.scatter(timesteps, MAEs, color='r')
    plt.plot(timesteps, MAEs, color='r')
    plt.scatter(timesteps, MAEs_s, color='b')
    plt.plot(timesteps, MAEs_s, color='b')

    plt.legend(['cartesian', 'spherical'])
    
    plt.savefig('figures/Al_average_force_length_MAE_sphere.png')
    plt.show()

if __name__ == "__main__":
    sizes = ['big', 'small', 'smaller']
    for s in sizes:
        plot_energies(s)
        plot_forces(s)
    
    # håll separata tränings- och evalueringsdataset
    # säg håll ett test

    # träna explicit för simulering, mata in säg 10 observerade värden?
    # när vi tränar, optimera för att mata in modellens egna prediktioner
    # och komma så nära det riktiga datat som möjligt
    # formulera om som ett tidsserie problem? autoregressiv modell?
    # börja med enkel KRR, jämför med detta? parallelt jämföra med MTP
    # också?
