#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import random

if __name__ == "__main__":
    print("ML example")
    figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    #plt.xscale('log')
    plt.grid(True)

    time = np.arange(0, 6, 0.2)
    sine = np.sin(time)
    noise = np.random.normal(0, 0.1, sine.shape)
    noisy_sine = sine + noise
    
    data = list(zip(time, noisy_sine))
    random.shuffle(data)
    tr_data = data[:15]
    ts_data = data[15:]

    plt.plot(time, sine, color='black')
    plt.scatter(time, noisy_sine)

    plt.legend(["y=sin(x)", "Complete data set"])
    plt.savefig('figures/ML/pure.png')

    # Noisy data and divided data
    figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    #plt.xscale('log')
    plt.grid(True)
    
    plt.plot(time, sine, color="black")
    plt.scatter(*zip(*tr_data), color='red')
    plt.scatter(*zip(*ts_data), color='green')

    plt.legend(["y=sin(x)", "Training data", "Testing data"])
    plt.savefig('figures/ML/divided_and_shuffled.png')

    # polynomial regression
    # get training data
    x,y = zip(*tr_data)
    
    model_line = np.arange(0,6,0.001)
    for i in range(1,10):
        figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
        plt.title("")
        plt.xlabel("")
        plt.ylabel("")
        #plt.xscale('log')
        plt.grid(True)
        plt.ylim([-1.5, 1.5])

        model = np.poly1d(np.polyfit(x,y,i)) 
        plt.plot(time, sine, color="black")
        plt.plot(model_line, model(model_line), linestyle='--')
        plt.scatter(*zip(*tr_data), color='red')
        plt.scatter(*zip(*ts_data), color='green')
        
        plt.legend(["y=sin(x)", "Poly. fit degree: " + str(i),"Training data", "Testing data"])
        plt.savefig('figures/ML/degree'+str(i)+'.png')
