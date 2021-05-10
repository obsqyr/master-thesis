#!/usr/bin/env python3
import md_qml
import os

if __name__ == "__main__":
    #index = [1, 10, 100, 1000, 10000]
    index = [1]
    for i in index:
        id = md_qml.run_md('MTP', i, 'Al', '06', 0, offset=2000)
