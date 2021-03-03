#!/usr/bin/env python3
"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.build import bulk
from asap3 import LennardJones
from asap3 import EMT
from asap3 import Trajectory
import properties as pr
import qml
import numpy as np
import calculators as calcs

def generate_qml_potential():
    alpha = np.loadtxt("machine.txt")

def run_md():
    # Set up a crystal
    atom = bulk('Al', 'fcc', a=4.0479, cubic=True)
    atoms = atom*(2,2,2)
    #size = 10
    #atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Al", size=(size, size, size), pbc=True)
    N = len(atoms.get_chemical_symbols())

    # Describe the interatomic interactions with the Effective Medium Theory
    #atoms.calc = EMT()
    # Describe the interatomic interactions with Lennard Jones potential
    #atoms.calc = LennardJones([18], [0.010323], [3.40], rCut = 6.625, modified = True)
    # Describe the interatomic interactions with QML
    #atoms.calc = calcs.zero_calculator()
    atoms.calc = calcs.KRR_calculator('Al', 1000)
    #generate_qml_potential()
    #x = qml_calc.qml_calculator()
    #print(x.get_potential_energy())
    #print(qml_calc.qml_calculator().get_potential_energy())
    
    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, temperature_K = 300)

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    #traj = Trajectory(atoms.get_chemical_symbols()[0] + '.traj', 'w', atoms)
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
    #dyn.attach(traj.write, interval=10)

    temperatures = []
    def printenergy(t=temperatures, a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        t.append(ekin / (1.5 * units.kB))
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

    # Now run the dynamics
    dyn.attach(printenergy, interval=10)
    print(printenergy())
    dyn.run(2000)

    size = 10
    return temperatures, N, atoms, size

if __name__ == "__main__":
    temperatures, N, atoms, size = run_md()

    # calculate specific heat
    spec_heat = pr.specific_heat(temperatures, N, atoms, size) / 1000 # convert to KJ/K*kg
    print ("Specific heat " + str(atoms.symbols) + ": %.4f [kJ/(K*kg)]" % (spec_heat))
