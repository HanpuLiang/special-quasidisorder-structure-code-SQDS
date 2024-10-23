'''
@author: Han-Pu Liang
@date: 2023/07/10
@description: This code find the best match of the target correlation function of the target supercell. 
'''

import numpy as np
import os, random

target_lro = 0.0
target_sro = 0.0
N_iter = 100000
expand_matrix = np.array([[5, 0, 0], [0, 5, 0], [0, 0, 2]])
exchange_atoms = ['Zn', 'Sn']
fixed_atoms = ['P']
order_file_name = 'ZnSnP2-symmetry-cell.vasp'

SRO_critic = 2  # 0: mean, 1: max, 2: exp(-(d-d0))

def get_cur_score(mpi_id):
    with open(f'cur_tcf_{mpi_id}.out', 'r') as obj:
        ct = obj.readlines()
    tcf_all = np.array([float(item) for item in ct[0].split()])
    return tcf_all

def perform_sqs(out_atat, out_positions_symbols, A_index, B_index, order_A_index, occupies, norm_weights, mpi_id):
    choose_A = random.randint(0, A_index.shape[0]-1)
    choose_B = random.randint(0, B_index.shape[0]-1)
    A_position_symbol = out_positions_symbols[A_index[choose_A]].split()
    B_position_symbol = out_positions_symbols[B_index[choose_B]].split()
    A_position_symbol[3], B_position_symbol[3] = B_position_symbol[3], A_position_symbol[3]
    out_positions_symbols[A_index[choose_A]] = ' '.join(A_position_symbol)
    out_positions_symbols[B_index[choose_B]] = ' '.join(B_position_symbol)
    occupies[A_index[choose_A]] ^= 1
    occupies[B_index[choose_B]] ^= 1
    cur_lro = 2*(1-np.sum(occupies[order_A_index])/len(order_A_index)) - 1
    A_index[choose_A], B_index[choose_B] = B_index[choose_B], A_index[choose_A]
    
    cur_position_lines = '\n'.join(out_positions_symbols)
    cur_atat = f'''{out_atat}
{cur_position_lines}
'''
    with open(f'cur_sqs_{mpi_id}.out', 'w') as obj:
        obj.write(cur_atat)
    os.system(f'corrdump -2=4.2 -l rndstr.in -noe -c -s cur_sqs_{mpi_id}.out 2>warning 1>cur_tcf_{mpi_id}.out')
    cur_tcf = get_cur_score(mpi_id)
    if SRO_critic == 0:
        cur_sro = np.mean(np.abs(cur_tcf))
    elif SRO_critic == 1:
        cur_sro = np.max(np.abs(cur_tcf))
    elif SRO_critic == 2:
        cur_sro = np.sum(np.abs(cur_tcf) * norm_weights)

    return out_positions_symbols, A_index, B_index, occupies, cur_lro, cur_sro, cur_tcf

def calc_score(target, current, weight):
    score = np.sum(np.abs(np.array(target) - np.array(current)) * weight)
    return score

def get_supercell_by_expand(cell, positions, expand_matrix):
    total_size = expand_matrix[0,0] * expand_matrix[1,1] * expand_matrix[2,2]
    expand_positions = np.zeros((positions.shape[0]*total_size, 3))
    count = 0
    for i in range(expand_matrix[0,0]):
        for j in range(expand_matrix[1,1]):
            for k in range(expand_matrix[2,2]):
                expand_positions[count*positions.shape[0]:(count+1)*positions.shape[0]] = positions + np.sum(cell, axis=0)*np.array([i, j, k])
                count = count + 1
    expand_cell = np.matmul(cell, expand_matrix)
    return expand_cell, expand_positions

def get_poscar_position(positions_symbol, A_index, B_index, other_index):
    A_positions = positions_symbol[A_index]
    B_positions = positions_symbol[B_index]
    all_positions = [A_positions, B_positions]
    for id in other_index:
        all_positions.append(positions_symbol[id])
    all_positions = list(np.concatenate(all_positions, axis=0))
    return '\n'.join(all_positions)
    


def convert_positions(cell, positions):
    cell = np.array([[float(item) for item in line.split()] for line in cell])
    positions = np.array([[float(item) for item in line.split()] for line in positions])
    new_positions = np.multiply(cell, positions) # this multiply may is not correct
    return new_positions 
    

def sqs_main():
    cur_mpi_id = os.getenv('SLURM_PROCID')
    
    # the desired correlation function
    score = 9999
    # the weight of each cluster, of course, an uniform weight is ok
    weight = np.array([0.6, 0.4]) ### <<< ----- change --- <<<
    expand_coe = expand_matrix[0,0] * expand_matrix[1,1] * expand_matrix[2,2]

    with open('clusters.out', 'r') as obj:
        ct_cluster = obj.readlines()
    cluster_dis = []
    for i in range(100):
        if 6*i+1 >= len(ct_cluster):
            break
        cluster_dis.append(float(ct_cluster[6*i+1].split()[0]))
    cluster_dis = np.array(cluster_dis)
    weights = np.exp(-(cluster_dis - np.min(cluster_dis)))
    norm_weights = weights / np.sum(weights)
    
    with open(order_file_name, 'r') as obj:
        order_poscar = obj.readlines()
    order_cell = np.array([[float(item) for item in line.split()] for line in order_poscar[2:5]])
    order_positions = np.array([[float(item) for item in line.split()[0:3]] for line in order_poscar[8:]])
    order_positions = np.matmul(order_positions, order_cell)
    atom_symbol = np.array([item for item in order_poscar[5].split()])
    atom_nums = np.array([float(item) for item in order_poscar[6].split()])
    atom_symbols = []
    for i, n in enumerate(atom_nums):
        for j in range(int(n)):
            atom_symbols.append(atom_symbol[i])

    super_cell, super_positions = get_supercell_by_expand(order_cell, order_positions, expand_matrix)
    super_symbols = np.array(atom_symbols * expand_coe)
    super_positions_frac = np.matmul(super_positions, np.linalg.inv(super_cell))
    super_positions_frac_str = [' '.join([f'{item:15.10f}' for item in super_positions_frac[i]] + [super_symbols[i]]) for i in range(len(super_positions))]

    A_index = np.where(super_symbols == exchange_atoms[0])[0]
    B_index = np.where(super_symbols == exchange_atoms[1])[0]
    order_A_index = A_index.copy()
    other_index = [np.where(super_symbols == fixed_atoms[i])[0] for i in range(len(fixed_atoms))]

    occupies = np.zeros(super_symbols.shape)
    occupies = occupies.astype(np.int16)
    occupies[B_index] = 1
    
    # get the former of poscar
    cell_vasp = '\n'.join([' '.join([f'{item:20.15f}' for item in line]) for line in super_cell])
    symbols_vasp = ' '.join([f'{item:10s}' for item in atom_symbol])
    number_vasp = ' '.join([f'{item*expand_coe:10.0f}' for item in atom_nums])
    out_poscar = f'''Supercell for SQDS
1.0
{cell_vasp}
{symbols_vasp}
{number_vasp}
Direct
'''
    out_atat = f'''{cell_vasp}
1 0 0
0 1 0
0 0 1'''

    out_poscar_positions = get_poscar_position(np.array(super_positions_frac_str), A_index, B_index, other_index)
    out_poscar_order = out_poscar + out_poscar_positions
    with open(f'expand_order_POSCAR.vasp', 'w') as obj:
        obj.write(out_poscar_order)

    # random structures
    print(' --> Start random structure')
    for i in range(100):
        print(f'     Random structure : {i}')
        cur_positions, cur_A_index, cur_B_index, cur_occupies, cur_lro, cur_sro, cur_tcf = perform_sqs(out_atat, super_positions_frac_str.copy(), A_index.copy(), B_index.copy(), order_A_index.copy(), occupies.copy(), norm_weights, cur_mpi_id)
        super_positions_frac_str = cur_positions
        A_index = cur_A_index
        B_index = cur_B_index
        occupies = cur_occupies

    for i in range(N_iter):
        cur_positions, cur_A_index, cur_B_index, cur_occupies, cur_lro, cur_sro, cur_tcf = perform_sqs(out_atat, super_positions_frac_str.copy(), A_index.copy(), B_index.copy(), order_A_index.copy(), occupies.copy(), norm_weights, cur_mpi_id)
        cur_score = calc_score([target_lro, target_sro], [cur_lro, cur_sro], weight)

        p = np.exp(-cur_score*150)

        if cur_score < score:
            score = cur_score
            super_positions_frac_str = cur_positions.copy()
            A_index = cur_A_index.copy()
            B_index = cur_B_index.copy()
            occupies = cur_occupies.copy()
            print(occupies)
            str_cur_lro_sro_tcf = ' '.join([f'{cur_lro:15.10f}, '] + [f'{cur_sro:15.10f}, '] + [f'{item:15.10f}' for item in cur_tcf])
            current_info = f'     Current lro, sro, tcf : {str_cur_lro_sro_tcf} ; Score : [ {cur_score:15.10f} ] \n'
            out_poscar_positions = get_poscar_position(np.array(super_positions_frac_str), A_index, B_index, other_index)
            cur_out_poscar = out_poscar + out_poscar_positions
            with open(f'best_POSCAR_{cur_mpi_id}.vasp', 'w') as obj:
                obj.write(cur_out_poscar)
            os.system(f'cp cur_sqs_{cur_mpi_id}.out best_sqs_{cur_mpi_id}.out')
            os.system(f'cp cur_tcf_{cur_mpi_id}.out best_tcf_{cur_mpi_id}.out')
            if score < 0.01:
                os.system(f'cp best_POSCAR_{cur_mpi_id}.vasp save-best-data/POSCAR-{score:.10f}-{cur_mpi_id}.vasp')
                os.system(f'cp best_sqs_{cur_mpi_id}.out save-best-data/sqs-{score:.10f}-{cur_mpi_id}.out')
                os.system(f'cp best_tcf_{cur_mpi_id}.out save-best-data/tcf-{score:.10f}-{cur_mpi_id}.out')
        elif random.random() <= p:
            current_info = f'     Metroplis rising !  Current lro, sro, tcf : {str_cur_lro_sro_tcf} ; Score : [ {cur_score:15.10f} ]'
            score = cur_score
            super_positions_frac_str = cur_positions
            A_index = cur_A_index
            B_index = cur_B_index
            occupies = cur_occupies
        print(current_info)
        with open(f'log-{cur_mpi_id}', 'a') as obj:
            obj.write(current_info)



if __name__ == '__main__':
    sqs_main() 
