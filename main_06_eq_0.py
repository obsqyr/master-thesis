#!/usr/bin/env python3
import md_qml
import time
import os

if __name__ == "__main__":
    if os.path.exists('property_calculations/eq_0/times.txt'):
        os.remove('property_calculations/eq_0/times.txt')
    f = open('property_calculations/eq_0/times.txt', 'w+')
    f.write('id\ttime[s]\n')

    index = [1, 10, 100, 1000, 10000]
    #index = [0,1,2,3,4,5,6,7,8,9]
    for i in index:
        #print(i)
        start_time = time.time()
        id = md_qml.run_md('MTP', i, 'Si', '06', 0, 'eq_0/', 0)
        #for simulating the same potential iteratively
        #id = md_qml.run_md('MTP', 100, 'Al', '06', 0, 'eq_0_iter_'+str(i)+'/')
        print("Ran in %s seconds" % (time.time() - start_time))
        f.write(id + '\t' + str(time.time() - start_time) + '\n')
    f.close()
    

