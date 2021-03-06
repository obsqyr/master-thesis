#!/usr/bin/env python3
import md_qml
import time
import os

if __name__ == "__main__":
    if os.path.exists('property_calculations/eq_2000/times.txt'):
        os.remove('property_calculations/eq_2000/times.txt')
    f = open('property_calculations/eq_2000/times.txt', 'w+')
    f.write('id\ttime[s]\n')

    index = [1, 10, 100, 1000, 10000]
    for i in index:
        start_time = time.time()
        id = md_qml.run_md('MTP', i, 'Al', '10', 2000, 'eq_2000/')
        print("Ran in %s seconds" % (time.time() - start_time))
        f.write(id + '\t' + str(time.time() - start_time) + '\n')
    f.close()

