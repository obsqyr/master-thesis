#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import os
import math

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


def plot_mtp_log():
    dirs = [x[0] for x in os.walk('test_results_log')]
    print(dirs)
    dirs.remove('test_results_log')
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
        
    for i, pot in enumerate(filenames_al):
        tuples_al = []
        for f in pot:
            if f[23:27].isdecimal():
                tuples_al.append((f, int(f[23:27])))
            elif f[23:26].isdecimal():
                tuples_al.append((f, int(f[23:26])))
            elif f[23:25].isdecimal():
                tuples_al.append((f, int(f[23:25])))
            elif f[23:24].isdecimal():
                tuples_al.append((f, int(f[23:24])))
        
            tuples_al = sorted(tuples_al, key=lambda x: x[1])
            filenames_al[i] = [i[0] for i in tuples_al]

    for i, pot in enumerate(filenames_si):
        tuples_si = []
        for f in pot:
            if f[23:27].isdecimal():
                tuples_si.append((f, int(f[23:27])))
            elif f[23:26].isdecimal():
                tuples_si.append((f, int(f[23:26])))
            elif f[23:25].isdecimal():
                tuples_si.append((f, int(f[23:25])))
            elif f[23:24].isdecimal():
                tuples_si.append((f, int(f[23:24])))
        
            tuples_si = sorted(tuples_si, key=lambda x: x[1])
            filenames_si[i] = [i[0] for i in tuples_si]

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
#    del MAEs_al[0]
    
    # plot data_al
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    log = np.logspace(0.1, 2.5, 50)
    log = [math.ceil(i) for i in log]
    #print("generating training .cfg files")
    timesteps = sorted(set(log))
    #timesteps = range(100,1000,100)

    # energy
    for pot in MAEs_al:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps[:10], MAE[:10])
        plt.plot(timesteps[:10], MAE[:10])

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Al_energy_per_atom_MAE_mtp_log.png')
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
        plt.scatter(timesteps[:10], MAE[:10])
        plt.plot(timesteps[:10], MAE[:10])

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Al_forces_MAE_mtp_log.png')
    #plt.show()

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
    print('MAEs_si', MAEs_si)
    #del MAEs_si[0]
    #del MAEs_al[1]
    
    # plot data_si_small
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    log = np.logspace(0.1, 2.5, 50)
    log = [math.ceil(i) for i in log]
    #print("generating training .cfg files")
    timesteps = sorted(set(log))
    
    # energy
    for pot in MAEs_si:
        MAE = [m[1] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps[:30], MAE[:30])
        plt.plot(timesteps[:30], MAE[:30])

    plt.legend(['06.mtp','10.mtp'])
    plt.savefig('figures/Si_energy_per_atom_MAE_mtp_log.png')
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
        plt.scatter(timesteps[:30], MAE[:30])
        plt.plot(timesteps[:30], MAE[:30])

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_forces_MAE_mtp_log.png')
    #plt.show()

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
    al_en = figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
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
    al_fo = figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
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
    si_en = figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
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
    si_fo = figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
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


def plot_mtp_cv():
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
        filenames_al_small.append([i for i in filename if i[19:21].isdecimal() and not i[19:22].isdecimal()]) 
    
    print(len(filenames_al_small))
    
    filenames_si_small = []
    for filename in filenames_si:
        filenames_si_small.append([i for i in filename if i[19:21].isdecimal() and not i[19:22].isdecimal()])

    filenames_al_small[1].append('test_results/06/Al_100.txt')
    filenames_al_small[2].append('test_results/10/Al_100.txt')
    filenames_si_small[1].append('test_results/06/Si_100.txt')
    filenames_si_small[2].append('test_results/10/Si_100.txt')
    print('filenames_al', filenames_al_small)
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

    MAEs_al_small = []
    for d in data_al_small:
        MAEs_al_small.append([extract_MAE(x) for x in d])
    #print(len(MAEs_al[2]))
    #MAEs_al.append(MAEs_al[0])
    del MAEs_al_small[0]
    #del MAEs_al[1]

    # Silicon
    data_si_small = []
    for filename in filenames_si_small:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si_small.append(d)

    MAEs_si_small = []
    for d in data_si_small:
        MAEs_si_small.append([extract_MAE(x) for x in d])
    #print('MAEs_si', MAEs_si)
    del MAEs_si_small[0]
    #del MAEs_al[1]

    # get cv results
    dirs = [x[0] for x in os.walk('cfg_cv/results')]
    dirs.remove('cfg_cv/results')
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

    filenames_al[0].append(filenames_al[0].pop(1))
    filenames_al[1].append(filenames_al[1].pop(1))
    #print(filenames_al)
    filenames_si[0].append(filenames_si[0].pop(1))
    filenames_si[1].append(filenames_si[1].pop(1))

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
    
    #print(MAEs_al)
    #print(len(MAEs_al[2]))
    #MAEs_al.append(MAEs_al[0])
    #del MAEs_al[0]
    #del MAEs_al[1]
    
    # plot data_al
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, 10-fold CV (Al, MTP)")
    plt.xlabel('timestep used as test data')
    plt.ylabel('MAE [eV]')
    timesteps = range(10,110,10)
    plt.xticks(timesteps)
    
    # energy
    #for pot in MAEs_al_small:
    MAE = [m[1] for m in MAEs_al_small[0]]
    print('timesteps', len(timesteps), 'MAE', len(MAE))
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    c = ['blue', 'orange']
    #for i, pot in enumerate(MAEs_al):
    MAE = [m[1] for m in MAEs_al[0]]
    print("Al, pot " + str(0) + " : Mean MAE for energy:", str(sum(MAE)/10))
    #print(len(MAE), len(timesteps))
    for m in MAE:
        plt.scatter(100, m, color=c[0])
        #plt.plot(timesteps, MAE)
        
    plt.legend(['06.mtp'])
    plt.savefig('figures/Al_energy_cv.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, 10-fold CV (Al, MTP)")
    plt.xlabel('timestep used as test data')
    plt.xticks(timesteps)
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    for i, pot in enumerate(MAEs_al):
        MAE = [m[2] for m in pot]
        print("Al, pot " + str(i) +" : Mean MAE for forces:", str(sum(MAE)/10))
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mpt'])
    plt.savefig('figures/Al_forces_cv.png')
    #plt.show()
    
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
    #print('MAEs_si', MAEs_si)
    #del MAEs_si[0]
    
    # plot data_si_small
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, 10-fold CV (Si, MTP)")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(10,110,10)
    plt.xticks(timesteps)

    # energy
    for i, pot in enumerate(MAEs_si):
        MAE = [m[1] for m in pot]
        print("Si, pot " + str(i) + " : Mean MAE for energy " + str(sum(MAE)/10))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp','10.mtp'])
    plt.savefig('figures/Si_energy_cv.png')
    #plt.show()

    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    plt.xticks(timesteps)
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    for i, pot in enumerate(MAEs_si):
        MAE = [m[2] for m in pot]
        print("Si, pot " + str(i) + " : Mean MAE for forces " + str(sum(MAE)/10))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)

    plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_forces_MAE_cv.png')
    #plt.show()
   

if __name__ == "__main__":
    #plot_mtp_closer_to_zero()
    plot_mtp_cv()
    #plot_mtp_log()

    
