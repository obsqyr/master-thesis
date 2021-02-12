#!/usr/bin/env python3
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("plotting")
    
    # total amount of timesteps, divided into test/training by
    # proportion 50/50
    timesteps = [10, 100, 200, 300, 400, 500]
    MAE = [0.053403660380666906, 2.9544048360212387, 1.273760723746573, 1.3170717980190503, 1.3331982309170984, 1.3232172627262855]

    plt.title("MAE against timesteps (divided into training/test by 50/50)")
    plt.xlabel('timesteps')
    plt.ylabel('MAE')
    plt.scatter(timesteps, MAE)
    plt.plot(timesteps, MAE)

    plt.savefig('figures/MAE_pos-pot.png')
    
