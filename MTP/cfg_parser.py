#!/usr/bin/env python3

# allowed elements
elements = ['Al', 'Si']

def generate_train_cfg(element, num_timesteps):
    if element not in elements:
        raise ValueError("Element " + element + " not allowed.")

    print("Generating train .cfg for " + element + " with first " + str(num_timesteps) + " timesteps")
    
    write_f = open('cfg_train/'+element+'_train_'+str(num_timesteps)+'.cfg', 'w+')

    num_cfgs = 0
    with open("cfg_out/"+element+"_relax.cfg") as f:
        for i, line in enumerate(f):
            stripped_line = line.strip()
            write_f.write(line)
            if line.strip() == "END_CFG":
                num_cfgs += 1
            if num_cfgs == num_timesteps:
                break
    write_f.close()
    

def generate_test_cfg(element, num_timesteps):
    if element not in elements:
        raise ValueError("Element " + element + " not allowed.")

    print("Generating test .cfg for " + element + " with last " + str(num_timesteps) + " timesteps")
    
    write_f = open('cfg_test/'+element+'_test_'+str(num_timesteps)+'.cfg', 'w+')

    num_cfgs = 0
    with open("cfg_out/"+element+"_relax.cfg") as f:
        for i, line in enumerate(f):
            stripped_line = line.strip()
            if line.strip() == "END_CFG":
                num_cfgs += 1
            # hard coded 10000, total amount of timesteps
            if num_cfgs >= (10000 - num_timesteps):
                write_f.write(line)
    # would like to remove first two lines in resulting file
    write_f.close()

    # very inefficient, rethink this at some point
    with open('cfg_test/'+element+'_test_'+str(num_timesteps)+'.cfg','r') as fin:
        data = fin.read().splitlines(True)
    with open('cfg_test/'+element+'_test_'+str(num_timesteps)+'.cfg', 'w') as fout:
        fout.writelines(data[2:])

if __name__ == "__main__":
    print("generating training .cfg files")
    for i in range(1000,10000,1000):
        generate_train_cfg('Al', i)
        generate_train_cfg('Si', i)
    print("generating testing .cfg files")
    generate_test_cfg('Al', 1000)
    generate_test_cfg('Si', 1000)

