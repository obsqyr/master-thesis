#!/usr/bin/env python3
import md_qml
import time
import os

if __name__ == "__main__":
    if os.path.exists('property_calculations/eq_2000/times.txt'):
        os.remove('property_calculations/eq_2000/times.txt')
    f = open('property_calculations/eq_2000/times.txt', 'w+')
    f.write('id\ttime[s]\n')

    index = [0,1,2,3,4,5,6,7,8,9]
    for i in index:
        #start_time = time.time()
        #id = md_qml.run_md('MTP', i, 'Al', '06', 2000, 'eq_2000/', 6000)
        # for training on the same potential iteratively
        id = md_qml.run_md('MTP', 100, 'Al', '06', 2000, 'eq_2000_iter_'+str(i)+'/', 2000)
        #print("Ran in %s seconds" % (time.time() - start_time))
        #f.write(id + '\t' + str(time.time() - start_time) + '\n')
    f.close()
