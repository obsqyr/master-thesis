#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os

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

if __name__ == "__main__":
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

    # håll separata tränings- och evalueringsdataset
    # säg håll ett test

    # träna explicit för simulering, mata in säg 10 observerade värden?
    # när vi tränar, optimera för att mata in modellens egna prediktioner
    # och komma så nära det riktiga datat som möjligt
    # formulera om som ett tidsserie problem? autoregressiv modell?
    # börja med enkel KRR, jämför med detta? parallelt jämföra med MTP
    # också?
