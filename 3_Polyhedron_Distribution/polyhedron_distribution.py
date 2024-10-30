'''
@author: Han-Pu Liang
@date: 2024/10/20
@affiliation: Beijing computational science research center; Eastern institute of technology, Ningbo
@email: hanpuliang@csrc.ac.cn
@description: This code calculates the polyhedra distribution in disordered alloys. 
                In the case of zinc-blende lattice ZnSnP2, the nearest-neighbor is tetrahedron, so the class number of tetrahedron is 5;
                the second-nearest-nerghbor polyhedron has 12 vertices, so the class number of polyhedron is 13.
                For the rock-salt lattice, the nearest-neighbor is octahedron.
'''

import numpy as np
from ase.io import read

def output_data(name, data):
    out_data = name + ' '.join([f'{i:6.2f}' for i in data])
    print(out_data)

''' calculate the polyhedra distribution of structure '''
def calc_distribution():
    mats = ['random-64atoms.vasp', 'random-512atoms.vasp']
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
        tetra_N_coe = np.sum(all_tetra[i,:]) / 16
        poly_N_coe = np.sum(all_poly_sec[i,:]) / 256
        output_data(mats[i]+' Tetra: ', all_tetra[i,:] / tetra_N_coe)
        output_data(mats[i]+' Poly-sec: ', all_poly_sec[i,:]/ poly_N_coe)

if __name__ == '__main__':
    calc_distribution()
