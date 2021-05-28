#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import os
import train_qml as tr
import vasp_parser as vp
import properties as pr

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
        'all' - MAEs for timesteps 10, 20, ..., 100, 200, ..., 9000
        'big' - MAEs for timesteps 1000, 2000, ... , 9000
        'small' - MAEs for timesteps 100, 200, ... , 900
        'smaller' - MAEs for timesteps 10, 20, ... , 90
        'smallest' - MAEs for timesteps 1, 2, ..., 9
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
        filenames_si = temp_si
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
    elif size == 'smallest':
        temp_al = []
        for filename in filenames_al:
            #print(filename)
            temp_al.append([i for i in filename if i[23:24].isdecimal() and not i[23:25].isdecimal()])
        filenames_al = temp_al
        temp_si = []
        for filename in filenames_si:
            #print(filename)
            temp_si.append([i for i in filename if i[23:24].isdecimal() and not i[23:25].isdecimal()])
        filenames_si = temp_si
    elif size == 'all':
        temp_al = []
        big = []
        small = []
        smaller = []
        smallest = []
        for filename in filenames_al:
            #print(filename)
            big.append([i for i in filename if i[23:27].isdecimal()])
            small.append([i for i in filename if i[23:26].isdecimal() and not i[23:27].isdecimal()])
            smaller.append([i for i in filename if i[23:25].isdecimal() and not i[23:26].isdecimal()])
            smallest.append([i for i in filename if i[23:24].isdecimal() and not i[23:25].isdecimal()])
        print('smallest', smallest)
        for i, element in enumerate(big):
            l = [smallest[i], smaller[i], small[i], big[i]]
            l = [item for sublist in l for item in sublist]
            temp_al.append(l)
        filenames_al = temp_al
        temp_si = []
        big = []
        small = []
        smaller = []
        smallest = []
        for filename in filenames_si:
            #print(filename)
            big.append([i for i in filename if i[23:27].isdecimal()])
            small.append([i for i in filename if i[23:26].isdecimal() and not i[23:27].isdecimal()])
            smaller.append([i for i in filename if i[23:25].isdecimal() and not i[23:26].isdecimal()])
            smallest.append([i for i in filename if i[23:24].isdecimal() and not i[23:25].isdecimal()])
        for i, element in enumerate(big):
            l = [smallest[i], smaller[i], small[i], big[i]]
            l = [item for sublist in l for item in sublist]
            temp_si.append(l)
        filenames_si = temp_si
    
    print('filenames_al', filenames_al)
    
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
    fig = figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    plt.xscale('log')
    plt.ylim([0,0.04])
    # this depends on size
    if size == 'big':
        timesteps = range(1000,10000,1000)
    elif size == 'small':
        timesteps = range(100, 1000, 100)
    elif size == 'smaller':
        timesteps = range(10, 100, 10)
    elif size == 'smallest':
        timesteps = range(1, 10, 1)
    elif size == 'all':
        t0 = list(range(1, 10, 1))
        t1 = list(range(10, 100, 10))
        t2 = list(range(100, 1000, 100))
        t3 = list(range(1000, 10000, 1000))
        timesteps = [t0, t1, t2, t3]
        # flatten list
        timesteps = [item for sublist in timesteps for item in sublist]
        
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
    #set_xscale('log')
    plt.savefig('figures/Al_energy_per_atom_MAE_all_' + size + '.png')
    plt.show()

    # plot Si energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Energy / atom MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV]')
    plt.ylim([0,0.16])
    plt.xscale('log')
    
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
    plt.title("Average force length MAE, Al")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV/ Å]')
    plt.xscale('log')
    plt.ylim([0,3.2])
    
    # this depends on size
    if size == 'big':
        timesteps = range(1000,10000,1000)
    elif size == 'small':
        timesteps = range(100, 1000, 100)
    elif size == 'smaller':
        timesteps = range(10, 100, 10)
    elif size == 'smallest':
        timesteps = range(1, 10, 1)
    elif size == 'all':
        t0 = list(range(1,10,1))
        t1 = list(range(10, 100, 10))
        t2 = list(range(100, 1000, 100))
        t3 = list(range(1000, 10000, 1000))
        timesteps = [t0, t1, t2, t3]
        # flatten list
        timesteps = [item for sublist in timesteps for item in sublist]
    
    # forces
    for pot in MAEs_al:
        MAE = [m[2] for m in pot]
        print(MAE)
        print(len(MAE), len(timesteps))
        plt.scatter(timesteps, MAE)
        plt.plot(timesteps, MAE)
    
    qml_MAEs_al = np.loadtxt('forces_MAEs/Al_' + size + '.txt')
    print(len(qml_MAEs_al))
    plt.scatter(timesteps, qml_MAEs_al)
    plt.plot(timesteps, qml_MAEs_al)

    plt.legend(['06.mtp', '10.mtp', 'KRR'])
    #plt.grid(True)
    plt.savefig('figures/Al_mean_forces_lengths_MAE_all_' + size + '.png')
    plt.show()

    # plot Si energy
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Average force length MAE, Si")
    plt.xlabel('timesteps')
    plt.ylabel('MAE [eV / Å]')
    plt.xscale('log')
    plt.ylim([0,5.7])
    
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

def plot_properties_convergence(element, eq, mtp, final=False, offset=0, variance=False):
    '''
    Plot time averages of properties to check whether or not
    they converge. 
    
    Params:
    eq (str): at which timestep point equilibrium is assumed to be
    achieved for MTP potential. 
    '''
    # extract eq. timestep
    eq_ts = int(eq[3:])
    
    # MSD
    figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
    plt.title("MSD convergence, " + element + ", MTP " + mtp +", eq. at " + str(eq_ts) + ", offset " + str(offset))
    plt.xlabel('timesteps')
    plt.ylabel('MSD [Å^2]')
    #plt.xscale('log')
    plt.grid(True)

    timesteps = range(100,10100, 100)
    # get DFT averages
    MSDs, Cvs, E_tots = pr.get_averaged_properties('properties_'+element+'_DFT_eq_'+ str(offset) +'.txt')
    MSDs_mtp = []
    Cvs_mtp = []
    E_tots_mtp = []
    
    if final:
        if element == 'Al':
            indeces = [1000, 10000]
            indeces_str = [str(i) for i in indeces]
        elif element == 'Si':
            indeces = [10, 100, 1000]
            indeces_str = [str(i) for i in indeces]
        legend = ['DFT']
        legend.extend(indeces_str)
    else:
        indeces = range(10, 100, 10)

    # get averages from indeces
    for i in indeces:
        MSD, Cv, E_tot = pr.get_averaged_properties(eq+'/properties_'+element+'_MTP_'+mtp+'_'+str(i)+'_'+eq+'_offset_'+str(offset)+'_ranfor_10000.txt', element)
        MSDs_mtp.append(MSD)
        Cvs_mtp.append(Cv)
        E_tots_mtp.append(E_tot)
        #print(eq+'/properties_'+element+'_MTP_'+mtp+'_'+str(i)+'_'+eq+'_ranfor_10000.txt')
        
    # get averages for mtp trained for 100 timesteps many times 
    # !! UPDATE THIS TO USE OFFSET !! 
    if variance:
        MSDs_100 = []
        Cvs_100 = []
        E_tots_100 = []
        if element == 'Al':
            j = '100'
        elif element == 'Si':
            j = '1000'
        for i in range(0,10,1):
            MSD, Cv, E_tot = pr.get_averaged_properties(eq+'_iter_'+str(i)+'/properties_'+element+'_MTP_'+mtp+'_'+j+'_eq_'+str(eq_ts)+'_offset_'+str(offset)+'_ranfor_10000.txt')
            MSDs_100.append(MSD)
            Cvs_100.append(Cv)
            E_tots_100.append(E_tot)

    # MSD
    # plot DFT data
    #print('DFT MSDs', len(MSDs), 'timesteps', len(timesteps[int(offset/100):]))
    #plt.scatter(timesteps, MSDs)
    plt.plot(timesteps[int(offset/100):], MSDs, color='black', linewidth=3)

    #timesteps = range(0, 9900, 100)
    # all MTP MSDs
    for i, MSDs in enumerate(MSDs_mtp):
        #print('MTP MSDs', len(MSDs), 'timesteps', len(timesteps[int(offset/100):]))
        #plt.scatter(timesteps, MSDs[:-1])
        if element == 'Al':
            if mtp == '10' and eq_ts == 0:
                #if i == 2 or i == 3 or i == 4:
                if offset != 0:
                    plt.plot(timesteps[int(offset/100):], MSDs)
                else:
                    plt.plot(timesteps[int(offset/100):], MSDs[:-1])
            else:
                # why MSDs[-1]?
                if offset != 0:
                    plt.plot(timesteps[int(offset/100):], MSDs)
                else:
                    plt.plot(timesteps[int(offset/100):], MSDs[:-1])
        elif element == 'Si':
            if mtp == '10' and eq_ts == 2000:
                if i != 1:
                    plt.plot(timesteps, MSDs[:-1])
            else:
                if offset != 0:
                    plt.plot(timesteps[int(offset/100):], MSDs)
                else:
                    plt.plot(timesteps[int(offset/100):], MSDs[:-1])

    if variance:
        for MSD in MSDs_100:
            #plt.scatter(timesteps, Cv[:-1])
            if element == 'Al':
                plt.plot(timesteps[int(offset/100):], MSD, color='orange', alpha=0.3)
            elif element == 'Si':
                if mtp == '10' and eq_ts == 2000:
                    plt.plot(timesteps, MSD[:-1], color='green', alpha=0.3)
                else:
                    plt.plot(timesteps, MSD[:-1], color='red', alpha=0.3)
    
    #plt.grid(True)
    if final:
        if element == 'Al':
            if mtp == '10' and eq_ts == 0:
                plt.legend(['DFT', '1', '10', '100', '1000', '10000'])
            else:
                plt.legend(legend)
        elif element == 'Si':
            if mtp == '10' and eq_ts == 2000:
                plt.legend(['DFT', '1', '100', '1000', '10000'])
            else:
                plt.legend(legend)
            
        plt.savefig('figures/convergence/'+element+'_MSD_convergence_MTP_'+mtp+'_'+eq+'_offset_'+str(offset)+'_final.png')
    else:
        plt.legend(['DFT', '10', '20', '30', '40', '50', '60', '70', '80', '90'])
        plt.savefig('figures/convergence/'+element+'_MSD_convergence_MTP_'+mtp+'_'+eq+'.png')
    #plt.show()
    
    # specific heat
    figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Specific heat capacity convergence, "+ element +", MTP " +mtp+ ", eq. at " + str(eq_ts) + ", offset " + str(offset))
    plt.xlabel('timesteps')
    plt.ylabel('Specific heat capacity [J/(K*Kg)]')
    if element == 'Al':
        print('remember the ylim')
        #plt.ylim([460, 600])
    elif element == 'Si':
        print('remember the ylim')
        #plt.ylim([460, 950])
    #plt.xscale('log')
    plt.grid(True)
    timesteps = range(100,10100, 100)

    #print(len(timesteps[int(offset/100):]), len(Cvs))
    #plt.scatter(timesteps, Cvs)
    plt.plot(timesteps[int(offset/100):], Cvs, color='black', linewidth=3)

    for Cv in Cvs_mtp:
        #plt.scatter(timesteps, Cv[:-1])
        if offset != 0:
            plt.plot(timesteps[int(offset/100):], Cv)
        else:
            plt.plot(timesteps[int(offset/100):], Cv[:-1])

    if variance:
        for Cv in Cvs_100:
            #plt.scatter(timesteps, Cv[:-1])
            if element == 'Al':
                plt.plot(timesteps[int(offset/100):], Cv, color='orange', alpha=0.3)
            elif element == 'Si':
                plt.plot(timesteps, Cv[:-1], color='red', alpha=0.3)
    
    if final:
        if element == 'Al':
            plt.legend(legend)
        elif element == 'Si':
            plt.legend(legend)

        plt.savefig('figures/convergence/'+element+'_Cv_convergence_MTP_'+mtp+'_'+eq+'_offset_'+str(offset)+'_final.png')
    else:
        plt.legend(['DFT', '10', '20', '30', '40', '50', '60', '70', '80', '90'])
        plt.savefig('figures/convergence/'+element+'_Cv_convergence_MTP_'+mtp+'_'+eq+'.png')
    #plt.show()
    
    # E_tot
    figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Total energy, " + element + ", MTP " + mtp +", eq. at " + str(eq_ts) + ', offset ' + str(offset))
    plt.xlabel('timesteps')
    plt.ylabel('total energy [eV/atom]')
    if element == 'Al':
        plt.ylim([-3.70, -3.645])
    elif element == 'Si':
        print('ylim,etot')
        #plt.ylim([-8, 4])
    #plt.xscale('log')
    plt.grid(True)    
    timesteps = range(100,10100, 100)
    
    #plt.scatter(timesteps, E_tots)
    plt.plot(timesteps[int(offset/100):], E_tots, color='black', linewidth=3)

    for E_tot in E_tots_mtp:
        #plt.scatter(timesteps, E_tot[:-1])
        if offset != 0:
            plt.plot(timesteps[int(offset/100):], E_tot)
        else:
            plt.plot(timesteps[int(offset/100):], E_tot[:-1])
    
    if variance:
        for E_tot in E_tots_100:
            #print('plotting E_tots_100', E_tot[:-1])
            if element == 'Al':
                plt.plot(timesteps[int(offset/100):], E_tot, color='orange', alpha=0.3)
            elif element == 'Si':
                plt.plot(timesteps, E_tot[:-1], color='red', alpha=0.3)
    
    if final:
        if element == 'Al':
            plt.legend(legend)
        elif element == 'Si':
            plt.legend(legend)

        plt.savefig('figures/convergence/'+element+'_Etot_convergence_MTP_'+mtp+'_'+eq+'_offset_'+str(offset)+'_final.png')
    else:
        plt.legend(['DFT', '10', '20', '30', '40', '50', '60', '70', '80', '90'])
        plt.savefig('figures/convergence/'+element+'_Etot_convergence_MTP_'+mtp+'_'+eq+'.png')
    #plt.show()
    plt.close('all')

def plot_test(element, eq, mtp, final=False):
    # extract eq. timestep
    eq_ts = int(eq[3:])
    
    # MSD
    figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("MSD convergence, " + element + ", " + mtp +" eq. at " + str(eq_ts))
    plt.xlabel('timesteps')
    plt.ylabel('MSD [Å^2]')
    #plt.xscale('log')
    timesteps = range(100,10100, 100)
    MSDs, Cvs, E_tots = pr.get_averaged_properties('properties_'+element+'_DFT_eq_'+ str(eq_ts) +'.txt')
    MSDs_mtp = []
    Cvs_mtp = []
    E_tots_mtp = []
    
    if final:
        if element == 'Al':
            indeces = [1, 10, 100, 1000, 10000]
        elif element == 'Si':
            indeces = [1, 10, 100, 1000, 10000]
    else:
        indeces = range(10, 100, 10)

    # get averages from indeces
    for i in indeces:
        MSD, Cv, E_tot = pr.get_averaged_properties(eq+'/properties_'+element+'_MTP_'+mtp+'_'+str(i)+'_'+eq+'_ranfor_10000_test.txt')
        MSDs_mtp.append(MSD)
        Cvs_mtp.append(Cv)
        E_tots_mtp.append(E_tot)
        #print(eq+'/properties_'+element+'_MTP_'+mtp+'_'+str(i)+'_'+eq+'_ranfor_10000.txt')
        
    # get averages for mtp trained for 100 timesteps many times 
    MSDs_100 = []
    Cvs_100 = []
    E_tots_100 = []
    if element == 'Al':
        j = '100'
    elif element == 'Si':
        j = '1000'
    for i in range(0,10,1):
        MSD, Cv, E_tot = pr.get_averaged_properties(eq+'_iter_'+str(i)+'/properties_'+element+'_MTP_'+mtp+'_'+j+'_eq_'+str(eq_ts)+'_ranfor_10000.txt')
        MSDs_100.append(MSD)
        Cvs_100.append(Cv)
        E_tots_100.append(E_tot)

    # MSD
    # plot DFT data
    print('MSDs', len(MSDs), 'timesteps', len(timesteps[int(eq_ts/100):]))
    #plt.scatter(timesteps, MSDs)
    plt.plot(timesteps[int(eq_ts/100):], MSDs, color='black', linewidth=3)

    timesteps = range(2000, 9900, 100)
    # all MTP MSDs
    for i, MSDs in enumerate(MSDs_mtp):
        #print('MSDs', len(MSDs), 'timesteps', len(timesteps))
        #plt.scatter(timesteps, MSDs[:-1])
        if element == 'Al':
            if mtp == '10' and eq_ts == 0:
                if i == 2 or i == 3 or i == 4:
                    plt.plot(timesteps, MSDs[:-1])
            else:
                plt.plot(timesteps, MSDs[:-1])
        elif element == 'Si':
            if mtp == '10' and eq_ts == 2000:
                if i != 1:
                    plt.plot(timesteps, MSDs[:-1])
            else:
                plt.plot(timesteps, MSDs[:-1])

        #plt.grid(True)
    if final:
        if element == 'Al':
            if mtp == '10' and eq_ts == 0:
                plt.legend(['DFT', '100', '1000', '10000'])
            else:
                plt.legend(['DFT', '100'])
        elif element == 'Si':
            if mtp == '10' and eq_ts == 2000:
                plt.legend(['DFT', '1', '100', '1000', '10000'])
            else:
                plt.legend(['DFT', '1', '10', '100', '1000', '10000'])

        plt.savefig('figures/convergence/'+element+'_MSD_test.png')
    else:
        plt.legend(['DFT', '10', '20', '30', '40', '50', '60', '70', '80', '90'])
        plt.savefig('figures/convergence/'+element+'_MSD_convergence_MTP_new_'+mtp+'_'+eq+'.png')
    plt.show()

if __name__ == "__main__":
    #plot_properties_convergence('Si', 'eq_2000', '06', True)
    #plot_test('Al', 'eq_2000', '06', True)

    mtps = ['06']
    eqs = ['eq_0', 'eq_2000']
    offsets = [0, 2000]
    
    for mtp in mtps:
        for eq in eqs:
            for offset in offsets:
                print('mtp: ', mtp, '. eq: ', eq, '. offset: ', offset) 
                plot_properties_convergence('Al', eq, mtp, True, offset, False)
                #plot_properties_convergence('Si', eq, mtp, True, offset)

    #plot_forces('all')
    #plot_energies('all')
    #sizes = ['big', 'small', 'smaller']
    #for s in sizes:
    #    plot_energies(s)
    #    plot_forces(s)

    #plot_energies('all')
    #plot_forces('all')
    # håll separata tränings- och evalueringsdataset
    # säg håll ett test

    # träna explicit för simulering, mata in säg 10 observerade värden?
    # när vi tränar, optimera för att mata in modellens egna prediktioner
    # och komma så nära det riktiga datat som möjligt
    # formulera om som ett tidsserie problem? autoregressiv modell?
    # börja med enkel KRR, jämför med detta? parallelt jämföra med MTP
    # också?
