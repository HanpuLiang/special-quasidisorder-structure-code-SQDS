import numpy as np
from ase.io import read
import ase.io.vasp
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm #  用正态分布与数据分布做比较 
import os, shutil


''' write the tetrahedra distribution of structure with CF '''
def func2():
    # mats = ['poscar-221.vasp', 'poscar-222.vasp', 'poscar-331.vasp', 'lso0-sro0-332-2.vasp', 'poscar-442.vasp']
    mats = ['lso0.75-sro0.75-442-1-old.vasp', 'poscar-0.5-0.5-442.vasp', 'lso0.25-sro0.25-442.vasp']
    # id, CF[0], tetra
    anion = 'P'
    cation = 'Zn'
    all_tetra = np.zeros((len(mats),5))
    all_poly_sec = np.zeros((len(mats),13))
    for i in range(len(mats)):
        poscar = read(f'{mats[i]}', format='vasp')
        all_symbols = poscar.symbols
        anion_id = all_symbols == anion
        anion_number = np.sum(anion_id)
        all_dis = poscar.get_all_distances(mic=True)
        # calculate the anion atom dis
        anion_dis = all_dis[anion_id]
        for Si in range(anion_number):
            anion_sort_id = np.argsort(anion_dis[Si,:])
            anion_sort_dis = anion_dis[Si,anion_sort_id]
            anion_neighb_symb = all_symbols[anion_sort_id][1:5]         # nearest neighbor
            anion_sec_neighb_symb = all_symbols[anion_sort_id][17:29]   # second neighbor
            cur_N_cation = np.sum(anion_neighb_symb==cation)
            cur_N_sec_poly = np.sum(anion_sec_neighb_symb==cation)
            all_tetra[i,cur_N_cation] = all_tetra[i,cur_N_cation] + 1
            all_poly_sec[i,cur_N_sec_poly] = all_poly_sec[i,cur_N_sec_poly] + 1
        print(all_tetra[i,:])
        print(all_poly_sec[i,:])
    out_data = '\n'.join([', '.join([f'{item:7.4f}' for item in line]) for line in all_tetra])
    with open('CF-tetra.csv', 'w') as obj:
        obj.write(out_data)

if __name__ == '__main__':
    # func1()
    func2()
    # func3()
    # func4()
