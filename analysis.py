#!/usr/bin/env python3
from read_settings import read_settings_file
import os
from matplotlib import pyplot
import numpy as np
import math
import glob
import properties as pr
import warnings

def find_eq_lc(fnames):
    """
    Parameters:
    fnames (list): A list of strings of file path/names

    Returns:
    tuple: returns a tuple of lattice constant, bulk modulus,
    filename and interpolated lattice constant.
    """
    LCa = 9999
    LCb = 9999
    LCc = 9999
    Etot = 9999
    N = ""
    E_list = []
    LCa_list = []
    LCb_list = []
    LCc_list = []
    V_list = []
    for name in fnames:
        f = open(name, "r+")
        lines = f.read().split("\n")
        E = float(lines[-1].split()[2])
        # The strucutre of header is always the same.
        l_a = float(lines[11].split()[7])
        l_b = float(lines[11].split()[8])
        l_c = float(lines[11].split()[9])
        E_list.append(E)
        LCa_list.append(l_a)
        LCb_list.append(l_b)
        LCc_list.append(l_c)
        V_list.append(float(lines[11].split()[10]))

        if E < Etot:
            Etot = E
            LCa = l_a
            LCb = l_b
            LCc = l_c
            N = name
        f.close()

    settings = read_settings_file()
    n = LCa_list[0] * LCb_list[0] * LCc_list[0] * settings['supercell_size']**3 / V_list[0]
    oLCa = LCa_list[settings['LC_steps']] # Original lattice constant a.
    s_list = [x / oLCa for x in LCa_list]
    oV = V_list[settings['LC_steps']]
    p = np.polyfit(s_list, E_list, 2)
    LCia = 0
    LCib = 0
    LCic = 0
    if p[0] <= 0:
        warnings.warn("Dynamically unstable in this range. " + str(fnames))
        B = 0
        LCi = 0
    else:
        i_scaling = -p[1]/(2*p[0]) # interpolated lattice scaling factor
        LCia = i_scaling*LCa
        LCib = i_scaling*LCb
        LCic = i_scaling*LCc
        E_interp = np.polyval(p, i_scaling)
        V_interp = oV * i_scaling**3
        q = np.polyfit(V_list, E_list, 2)
        B = V_interp*(2*q[0])*160.2 # conversion from ev/Å^3 to GigaPa

    return LCa, LCb, LCc, B, N, LCia, LCib, LCic

def sort_properties_files():
    """
    Paramters:
    None

    Returns:
    tuple: returns a tuple of lists of lattice constants, bulk moduli,
    interpolated lattice constants and file paths.
    """
    settings = read_settings_file()
    filenames = sorted(glob.glob("property_calculations/properties_*"))
    steps = 1 + 2*settings['LC_steps']
    LC_list = []
    BulkM_list = []
    N_list = []
    LCi_list = []

    for i in range(0,round(len(filenames)/steps)):
        LCa, LCb, LCc, BulkM, N, LCia, LCib, LCic = find_eq_lc(filenames[steps*i:steps*(i+1)])
        LC_list.append([LCa, LCb, LCc])
        BulkM_list.append(BulkM)
        N_list.append(N)
        LCi_list.append([LCia, LCib, LCic])
        """
        for fname in filenames[steps*i:steps*i+steps]:
            if fname != N:
                os.remove(fname)
        """
    return LC_list, LCi_list, BulkM_list, N_list

def extract():
    """
    Paramters:
    None

    Returns:
    None
    """
    file = open("property_calculations/collected_data.txt", "w+")
    settings = read_settings_file()
    d = settings['decimals']

    def lj(str, k = d):
        return " "+str.ljust(k+10)

    file.write(lj("Material ID")+lj("Material")+lj("Cohesive energy")+lj("MSD")+lj("Self_diff")+lj("Specific heat"))

    if settings['vol_relax']:
        file.write(lj("Lattice const a")+lj("Lattice const b")+lj("Lattice const c"))
        file.write(lj("Interp LC a")+lj("Interp LC b")+lj("Interp LC c"))
        file.write(lj("Bulk modulus"))

    file.write(lj("Debye",2)+lj("Lindemann"))
    file.write("\n")

    file.write(lj(" ")+lj(" ")+lj("eV/atom")+lj("Å^2")+lj("mm^2/s")+lj("J/(K*Kg)"))

    if settings['vol_relax']:
        file.write(lj("Å")+lj("Å")+lj("Å")+lj("Å")+lj("Å")+lj("Å")+lj("Pa"))

    file.write(lj("K",2)+lj("1"))
    file.write("\n")

    file.close()
    N_list = glob.glob("property_calculations/properties_*")
    if settings['vol_relax']:
        LC_list, LCi_list, BulkM_list, N_list = sort_properties_files()

    for i, filename in enumerate(sorted(N_list)):
        f = open(filename, "r")
        lines = f.read().split("\n")
        f.close()
        if lines[-4] == 'Time averages:':
            matID = lines[0].split(":")[1]
            mat = lines[2].split()[1]
            Ecoh = lines[-1].split()[0]
            msd = lines[-1].split()[4]
            selfd = lines[-1].split()[5]
            Cv = lines[-1].split()[7]
            file = open("property_calculations/collected_data.txt", "a+")
            file.write(lj(matID)+lj(mat)+lj(Ecoh)+lj(msd)+lj(selfd)+lj(Cv))
            if settings['vol_relax']:
                LC = LC_list[i]
                LCi = LCi_list[i]
                BulkM = BulkM_list[i]
                file.write(pr.ss(LC[0],d+4)+pr.ss(LC[1],d+4)+pr.ss(LC[2],d+4))
                file.write(pr.ss(LCi[0],d+4)+pr.ss(LCi[1],d+4)+pr.ss(LCi[2],d+4))
                file.write(pr.ss(BulkM,d+4))
            if len(lines[-1].split()) > 8:
                debye = lines[-1].split()[8]
                linde = lines[-1].split()[9]
                file.write(pr.ss(debye,d+4)+pr.ss(linde,d+4))
            file.write("\n")
            file.close()
    return

def plot_properties():
    msd = []
    selfd = []
    spec_h = []
    latt_c = []
    inter_latt_c = []
    bulk_m = []
    coh_en = []
    debye = []
    linde = []

    f = open("property_calculations/collected_data.txt", "r")

    lines = f.readlines()[2:]
    for x in lines:
        coh_en.append(float(x.split()[2]))
        msd.append(float(x.split()[3]))
        selfd.append(float(x.split()[4]))
        spec_h.append(float(x.split()[5]))
        latt_c.append(float(x.split()[6]))
        inter_latt_c.append(float(x.split()[9]))
        bulk_m.append(float(x.split()[12]))
        if len(x.split()) > 13:
            debye.append(float(x.split()[13]))
            linde.append(float(x.split()[14]))

    bulk_m_filtered = []
    latt_c_filtered = []
    for i, value in enumerate(bulk_m):
        if value != 0 and abs(value) < 15000 and value > 0:
            bulk_m_filtered.append(value)
            latt_c_filtered.append(latt_c[i])

    f.close()
    #Plotting mean square displacment vs self diffusion const
    # in figure 1
    pyplot.figure(1)
    pyplot.scatter(msd,selfd)
    #Labeling the axes with names from properties.py
    pyplot.xlabel("Mean square displacement [Å^2]")
    pyplot.ylabel("Self diffusion [Å^2/fs]")

    pyplot.savefig("figures/MSD-SD.png")

    pyplot.figure(2)
    pyplot.scatter(latt_c_filtered,bulk_m_filtered)
    pyplot.xlabel("Lattice constant [Å]")
    pyplot.ylabel("Bulk modulus [GPa]")
    pyplot.savefig("figures/LC-BM.png")

    pyplot.figure(3)
    pyplot.scatter(latt_c, coh_en)
    pyplot.xlabel("Lattice constant [Å]")
    pyplot.ylabel("Cohesive energy [eV/atom]")
    pyplot.savefig("figures/LC-Ecoh.png")

    pyplot.figure(4)
    pyplot.scatter(latt_c, inter_latt_c)
    pyplot.xlabel("Lattice constant [Å]")
    pyplot.ylabel("Interpolated lattice constant [Å]")
    pyplot.savefig("figures/LC-inter_LC.png")

    pyplot.figure(5)
    pyplot.scatter(latt_c, spec_h)
    pyplot.xlabel("Lattice constant [Å]")
    pyplot.ylabel("Specific heat [J/(K*KG)]")
    pyplot.savefig("figures/LC-inter_LC.png")

    pyplot.show()

    return

if __name__ == "__main__":
    pr.clean_property_calculations()
    extract()
    #plot_properties()
