'''
@author: Han-Pu Liang
@date: 2023/07/10
@description: This code find the best match of the target correlation function of the target supercell. 
'''

import numpy as np
import os, random

def get_cur_score():
    with open('cur_tcf.out', 'r') as obj:
        ct = obj.readlines()
    tcf_all = np.array([float(item) for item in ct[0].split()])
    #tcf_id = np.array([1, 3])
    #tcf = tcf_all[tcf_id]
    return tcf_all

def perform_sqs(positions, A_site_id, B_site_id, cell, factor, element, N_atom):
    chose_A = random.choice(A_site_id)
    chose_B = random.choice(B_site_id)
    positions[chose_A], positions[chose_B] = positions[chose_B], positions[chose_A]
    
    # print(len(element))
    out_poscar = [f'{positions[i][:-1]} {element[i]} \n' for i in range(np.sum(N_atom))]
    out_poscar = ' '.join(cell+factor+out_poscar)
    with open('cur_sqs.out', 'w') as obj:
        obj.write(out_poscar)
    
    os.system('corrdump -2=4.2 -s=cur_sqs.out -l=rndstr.in -noe -c > cur_tcf.out')
    cur_tcf = get_cur_score()
    # more ctprint(cur_tcf)
    
    return positions, cur_tcf

def calc_score(cur_tcf, tcf, weight):
    score = np.sum(np.abs(cur_tcf - tcf)*np.array(weight))
    return score


def sqs_main():
    # the desired correlation function
    tcf = [-0.8503, 0.8503, -0.7723, 0.7046, 0.7115]  ### <<< ----- change --- <<<
    # the weight of each cluster, of course, an uniform weight is ok
    weight = [0.3, 0.3, 0.2, 0.1, 0.1] ### <<< ----- change --- <<<
    # file name of supercell
    file_name = 'MCS-443.vasp' ### <<< ----- change --- <<<
    
    
    
    with open(file_name, 'r') as obj:
        poscar = obj.readlines()
    cell = poscar[2:5]
    factor = ['1 0 0\n0 1 0\n0 0 1\n']
    element = [item for item in poscar[5].split()]
    N_atom = np.array([int(float(item)) for item in poscar[6].split()])
    N_total = np.sum(N_atom)
    positions = poscar[8:]
    
    # the index of changable atoms
    A_site_id = np.array(list(range(N_atom[0] + N_atom[1])))  ### <<< ----- change --- <<<
    
    
    B_site_id = A_site_id
    all_element = [element[0]]*N_atom[0] + [element[1]]*N_atom[1] + [element[2]]*N_atom[2]
    
    print(' --> Read structure sucessfully!')
    
    N_iter = 3000
    score = 9999
    
    
    
    
    
    if not os.path.exists('save-best-data'):
        print(' --> Create directory save-best-data/')
        os.makedirs('save-best-data')
    
    print(' --> Start random structure')
    for i in range(100):
        print(f'     Random structure : {i}')
        cur_positions, cur_tcf = perform_sqs(positions.copy(), A_site_id, B_site_id, cell, factor, all_element, N_atom)
        positions = cur_positions
    
    print(' --> Ready to perform Specical Quasi-disordered Structure')
    
    for i in range(N_iter):
        cur_positions, cur_tcf = perform_sqs(positions.copy(), A_site_id, B_site_id, cell, factor, all_element, N_atom)
        cur_score = calc_score(cur_tcf, tcf, weight)
        p = np.exp(-cur_score*100)
        
        if cur_score < score:
            score = cur_score
            positions = cur_positions
            str_cur_tcf = ' '.join([f'{item:10.6f}' for item in cur_tcf])
            print(f'     Current tcf : {str_cur_tcf} ; Score : [ {cur_score:10.6f} ]')
            out_poscar = ' '.join(poscar[0:8] + positions)
            with open('best_POSCAR.vasp', 'w') as obj:
                obj.write(out_poscar)
            os.system(f'cp cur_sqs.out best_sqs.out')
            os.system(f'cp cur_tcf.out best_tcf.out')
            os.system(f'cp best_POSCAR.vasp save-best-data/POSCAR-{score:.6f}.vasp')
            os.system(f'cp best_tcf.out save-best-data/tcf-{score:.6f}.out')
        elif random.random() < p:
            print(f'     Metroplis rising ! Current tcf : {str_cur_tcf} ; Score : [ {cur_score:10.6f} ]')
            score = cur_score
            positions = cur_positions

if __name__ == '__main__':
    sqs_main() 
