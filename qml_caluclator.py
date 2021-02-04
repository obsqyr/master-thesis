#!/usr/bin/env python3

from ase.calculators.calculator import Calculator

class qml_calculator(Calculator):
    def printname(self):
        print('qml')

if __name__ == "__main__":
    x = qml_calculator()
    x.printname()
