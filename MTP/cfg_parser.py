#!/usr/bin/env python3

from ase.io.cfg import *
from ase.io.vasp import *

# allowed elements
elements = ['Al', 'Si']

def atoms_to_cfg(atoms, fn):
    '''
    Function that writes atoms to .cfg format appropriate for MTP
    
    Restrictions:
    - Can only manage monoatomic materials
    '''
    print("converting atoms to mtp cfg")
    write_cfg('test.cfg', atoms)
    #write_vasp(fn, atoms)
    
    # create output string
    s = ""
    # write header information
    s += 'BEGIN_CFG\n' + ' Size\n' + '    ' + str(len(atoms.get_chemical_symbols())) + '\n Supercell\n'
    
    # write unit cell
    #print(atoms.get_cell())
    #if len(atoms.get_cell()) == 3:
    # orthorhombic cell
    for length in atoms.get_cell():
        space1 = '         '
        space2 = '      '
        space3 = '           '
        s += space1 + '{:.6f}'.format(length[0]) + space2 + '{:.6f}'.format(length[1]) + space2 + '{:.6f}'.format(length[2]) + '\n' 
        
    # write atom positions and forces
    s += ' AtomData:' + '  ' + 'id ' + 'type' + space1[:-2] + 'cartes_x' + space1[:-3] + 'cartes_y' + space1[:-3] + 'cartes_z' + space3 + 'fx' + space3[:-1] + 'fy' + space3[:-1] + 'fz\n'  
    for i, atom in enumerate(atoms):
        #print(atom.position)
        # currently only support for homogenous atoms
        pos = atom.position
        if i+1 < 10:
            s += '             ' + str(i+1)
        else:
            s += '            ' + str(i+1)
        s += '    ' + '0' + '       ' + '{:.6f}'.format(pos[0]) + '      ' + '{:.6f}'.format(pos[1]) + '      '+ '{:.6f}'.format(pos[0])+ '     ' +  '{:.6f}'.format(0) + '    ' + '{:.6f}'.format(0) + '    '+ '{:.6f}'.format(0)+ '\n'
            
    # write energy
    s += ' Energy\n' + '        0\n' 
    # write stress
    s += ' PlusStress:'+'  '+ 'xx'+ '          ' + 'yy' +'          '+ 'zz ' + '         '+ 'yz'+ '          ' + 'xz'+ '          '+ 'xy' +'\n\n'
    # write feature
    s += ' Feature\n'
    s += 'END_CFG\n'

    with open(fn, 'w+') as f:
        f.write(s)
   
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
    for i in range(10,100,10):
        generate_train_cfg('Al', i)
        generate_train_cfg('Si', i)
    #generate_train_cfg('Al', 1)
    #print("generating testing .cfg files")
    #generate_test_cfg('Al', 1000)
    #generate_test_cfg('Si', 1000)

