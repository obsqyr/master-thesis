#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.lines as mlines
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
    
    # plot data_si
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/atom]')
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

    line_06 = mlines.Line2D([], [], color='C0', marker='o', label='06.mtp')
    line_10 = mlines.Line2D([], [], color='C1', marker='o', label='10.mtp')
     
    plt.legend(handles=[line_06, line_10])
    #plt.legend(['06.mtp','10.mtp'])
    plt.grid(True)
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

    line_06 = mlines.Line2D([], [], color='C0', marker='o', label='06.mtp')
    line_10 = mlines.Line2D([], [], color='C1', marker='o', label='10.mtp')
     
    plt.legend(handles=[line_06, line_10])

    plt.grid(True)
    #plt.legend(['06.mtp', '10.mtp'])
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
    #print('filenames_al', filenames_al_small)
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
    dirs.remove('cfg_cv/results/06')
    dirs.remove('cfg_cv/results/10')
    dirs = sorted(dirs)
    dirs.remove('cfg_cv/results/06/100')
    dirs.append('cfg_cv/results/06/100')
    
    #print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)
    
    for f in filenames:
        x = f.pop(1)
        f.append(x)
        
    #print(filenames[1])
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
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    timesteps = range(10,110,10)
    plt.xticks(timesteps)
    
    # energy
    #for pot in MAEs_al_small:
    print(MAEs_al_small[0])
    MAE = [m[1] for m in MAEs_al_small[0]]
    print('timesteps', len(timesteps), 'MAE', len(MAE))
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    c = ['blue', 'orange']
    for i, pot in enumerate(MAEs_al):
        i += 1
        MAE = [m[1] for m in pot]
        print("Al, pot " + str(i) + " : Mean MAE for energy:", str(sum(MAE)/10))
        #print(len(MAE), len(timesteps))
        for m in MAE:
            j = i*10
            plt.scatter(j, m, color=c[1])
            #plt.plot(timesteps, MAE)
        
    plt.legend(['06.mtp'])
    plt.savefig('figures/Al_energy_cv.png')
    #plt.show()
    
    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE, 10-fold CV (Al, MTP)")
    plt.xlabel('timestep')
    plt.xticks(timesteps)
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    MAE = [m[2] for m in MAEs_al_small[0]]
    print('timesteps', len(timesteps), 'MAE', len(MAE))
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    for i, pot in enumerate(MAEs_al):
        i += 1
        MAE = [m[2] for m in pot]
        print("Al, pot " + str(i) +" : Mean MAE for forces:", str(sum(MAE)/10))
        for m in MAE:
            j = i*10
            plt.scatter(j, m, color=c[1])

    plt.legend(['06.mtp'])
    plt.savefig('figures/Al_forces_cv.png')
    #plt.show()

    '''
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
    '''
def plot_learning_curves_log():
    dirs_te = [x[0] for x in os.walk('test_results_log')]
    dirs_tr = [x[0] for x in os.walk('train_results_log')]
    #print(dirs)
    dirs_te.remove('test_results_log')
    dirs_tr.remove('train_results_log')
    #print(dirs)

    filenames_te = []
    for dir in dirs_te:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames_te.append(f)

    filenames_tr = []
    for dir in dirs_tr:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames_tr.append(f)

    #print(filenames_tr)

    # this could be automated
    #filenames_al = [filenames[0][0:9], filenames[1][0:9]]
    #filenames_si = [filenames[0][9:18], filenames[1][9:18]]
    #filenames_al = []
    filenames_si_te = []
    filenames_si_tr = []
    for filename in filenames_te:
        #filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si_te.append([f for f in filename if 'Si' in f])

    for filename in filenames_tr:
        #filenames_al.append([f for f in filename if 'Al' in f])
        filenames_si_tr.append([f for f in filename if 'Si' in f])
    
    '''
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
    '''
    
    for i, pot in enumerate(filenames_si_te):
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
            filenames_si_te[i] = [i[0] for i in tuples_si]

    #print(filenames_si_te)
            
    '''
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
    '''
    # Silicon
    data_si_te = []
    for filename in filenames_si_te:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si_te.append(d)
    
    data_si_tr = []
    for filename in filenames_si_tr:
        d = []
        for f in filename:
            with open(f, 'r') as f:
                d.append(f.read())
        data_si_tr.append(d)

    MAEs_si_te = []
    for d in data_si_te:
        MAEs_si_te.append([extract_MAE(x) for x in d])
    MAEs_si_tr = []
    for d in data_si_tr:
        MAEs_si_tr.append([extract_MAE(x) for x in d])
    
    #print('MAEs_si_te', MAEs_si_tr)
    #del MAEs_si[0]
    #del MAEs_al[1]
    
    # plot data_si
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE (Si, 06.mtp)")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/atom]')
    log = np.logspace(0.1, 2.5, 50)
    log = [math.ceil(i) for i in log]
    #print("generating training .cfg files")
    timesteps = sorted(set(log))
    
    # energy
    #for pot in MAEs_si:
    MAE_te = [m[1] for m in MAEs_si_te[0]]
    MAE_tr = [m[1] for m in MAEs_si_tr[0]]
    #print(MAE)
    #print(len(MAE), len(timesteps))
    plt.scatter(timesteps[:30], MAE_te[:30])
    plt.plot(timesteps[:30], MAE_te[:30])
    plt.scatter(timesteps[:30], MAE_tr[:30])
    plt.plot(timesteps[:30], MAE_tr[:30])

    line_te = mlines.Line2D([], [], color='C0', marker='o', label='test')
    line_tr = mlines.Line2D([], [], color='C1', marker='o', label='train')
     
    plt.legend(handles=[line_te, line_tr])
    #plt.legend(['06.mtp','10.mtp'])
    plt.grid(True)
    plt.savefig('figures/Si_energy_per_atom_learning_curves.png')
    #plt.show()
    
    # forces
    plt.clf()
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Forces MAE (Si, 06.mtp)")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/Å]')
    
    #MAE = []
    #MAE = [m[2] for m in MAEs]
    #print(MAE)
    #for pot in MAEs_si:
    MAE_te = [m[2] for m in MAEs_si_te[0]]
    MAE_tr = [m[2] for m in MAEs_si_tr[0]]
    #print(MAE)
    #print(len(MAE), len(timesteps))
    plt.scatter(timesteps[:30], MAE_te[:30])
    plt.plot(timesteps[:30], MAE_te[:30])
    plt.scatter(timesteps[:30], MAE_tr[:30])
    plt.plot(timesteps[:30], MAE_tr[:30])
     
    plt.legend(handles=[line_te, line_tr])

    plt.grid(True)
    #plt.legend(['06.mtp', '10.mtp'])
    plt.savefig('figures/Si_forces_MAE_learning_curves.png')
    #plt.show()

def plot_mtp_training_and_validation_errors():
    # read in iter files
    dirs = [x[0] for x in os.walk('training_output_iter')]
    #dirs.remove('test_results')
    print(dirs)

    filenames = []
    for dir in dirs:
        f = [dir+'/'+f for f in os.listdir(dir)]
        f = sorted(f)
        filenames.append(f)
        
    print(len(filenames))

    MAEs_si = []
    MAEs_al = []
    for filename in filenames[0]:
        print(filename)
        #print(filename[21:23])
        f = open(filename)
        i = f.read()
        if filename[21:23] == 'Si':
            MAEs_si.append(extract_MAE(i))
        elif filename[21:23] == 'Al':
            MAEs_al.append(extract_MAE(i))
        else:
            continue

    print(len(MAEs_si))
    print(len(MAEs_al))

    # ALUMINIUM
    figure(num=None, figsize=(8, 4.4), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy per atom average absolute difference, Al MTP 06")
    plt.xlabel('timesteps')
    plt.ylabel('Error [eV]')
    timesteps = [1, 10, 100, 1000, 8000]
    tr_er = [3.44613e-13, 2.82733e-05, 0.000246437, 0.000578199, 0.000501822]
    val_er = [0.208856, 0.0107632, 0.000615748, 0.000662907, 0.000595121]
    tr_er_2 = [1.83409e-13, 2.58173e-05, 0.000436566, 0.000517493, 0.00051914]
    val_er_2 = [0.0616174, 0.00512976, 0.000632198, 0.000577663, 0.000565452]
    plt.xscale('log')
    plt.yscale('log')
    
    tr_ers = []
    val_ers = []
    for m in MAEs_al:
        #plt.scatter(timesteps, m[1], color='blue')
        #plt.scatter(timesteps, m[1], color='blue')
        #plt.scatter(100, m[1], color='blue', alpha=0.3)
        #plt.scatter(100, m[6], color='orange', alpha=0.3)
        tr_ers.append(m[1])
        val_ers.append(m[6])

    plt.ylim([0.0001, 0.01])
    tr_std = np.std(tr_ers)
    val_std = np.std(val_ers)

    plt.scatter(timesteps, tr_er, color='blue')
    plt.scatter(timesteps, val_er, color='orange')
    plt.plot(timesteps, tr_er, tr_std, color='blue')
    plt.plot(timesteps, val_er, val_std, color='orange')
    plt.errorbar(100, tr_er[2], tr_std, capsize=3)
    plt.errorbar(100, val_er[2], val_std, capsize=3)

    #plt.scatter(timesteps, tr_er_2)
    #plt.scatter(timesteps, val_er_2)
    #plt.plot(timesteps, tr_er_2)
    #plt.plot(timesteps, val_er_2)
        
    plt.legend(['training', 'validation'])
    plt.savefig('figures/Al_errors.png')

    # SILICON
    figure(num=None, figsize=(8, 4.4), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy per atom average absolute difference, Si MTP 06")
    plt.xlabel('timesteps')
    plt.ylabel('Error [eV]')
    timesteps = [1, 10, 100, 1000, 8000]
    tr_er = [4.07674e-13, 2.83252e-06, 0.000155095, 0.000143081, 0.000120598]
    val_er = [0.059846, 0.00117006, 0.000269087, 0.000138267, 0.00010319]
    plt.xscale('log')
    plt.yscale('log')

    tr_ers = []
    val_ers = []
    for m in MAEs_si:
        #plt.scatter(timesteps, m[1], color='blue')
        #plt.scatter(timesteps, m[1], color='blue')
        #plt.scatter(100, m[1], color='blue', alpha=0.3)
        #plt.scatter(100, m[6], color='orange', alpha=0.3)
        tr_ers.append(m[1])
        val_ers.append(m[6])

    plt.ylim([0.00001, 0.01])
    tr_std = np.std(tr_ers)
    val_std = np.std(val_ers)

    plt.scatter(timesteps, tr_er, color='blue')
    plt.scatter(timesteps, val_er, color='orange')
    plt.plot(timesteps, tr_er, tr_std, color='blue')
    plt.plot(timesteps, val_er, val_std, color='orange')
    plt.errorbar(100, tr_er[2], tr_std, capsize=3)
    plt.errorbar(100, val_er[2], val_std, capsize=3)

    plt.legend(['training', 'validation'])
    plt.savefig('figures/Si_errors.png')


if __name__ == "__main__":
    #plot_mtp_closer_to_zero()
    #plot_mtp_cv()
    #plot_mtp_log()
    #plot_learning_curves_log()
    plot_mtp_training_and_validation_errors()
    
