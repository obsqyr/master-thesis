#!/usr/bin/env python3
import md_qml
import time
import os

if __name__ == "__main__":
    if os.path.exists('property_calculations/eq_0/times.txt'):
        os.remove('property_calculations/eq_0/times.txt')
    f = open('property_calculations/eq_0/times.txt', 'w+')
    f.write('id\ttime[s]\n')

    for i in range(10,110,10):
        start_time = time.time()
        id = md_qml.run_md('MTP', i, 'Al', '06', 0, 'eq_0/')
        print("Ran in %s seconds" % (time.time() - start_time))
        f.write(id + '\t' + str(time.time() - start_time) + '\n')
    f.close()
