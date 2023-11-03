'''
@author: Han-Pu Liang
@date: 2023/06/05
@description: This sqs2poscar code has no atom number limit. you can enter: `python sqs2poscar-HPLiang.py [file]` to convert the structure file of ATAT to POSCAR of VASP
'''

import sys
import numpy as np

def count_TRUE(array):
    num=0
    for i in range(len(array)):
        if array[i] == True: num+=1
    return num 

file_name = sys.argv[1]
with open(file_name, 'r') as obj:
    ct = obj.readlines()

cell = np.array([[float(item) for item in line.split()] for line in ct[0:3]])
lat = np.array([[float(item) for item in line.split()] for line in ct[3:6]])
pos = np.array([[float(item) for item in line.split()[0:3]] for line in ct[6:]])
atoms = np.array([line.split()[3] for line in ct[6:]]) 

decell = np.matmul(lat, cell)
inv_decell = np.linalg.inv(decell)
pos_frac = np.matmul( np.matmul(pos, cell), inv_decell )

atom_type = np.unique(atoms)
print(atom_type)
atom_num = []
pos_sort = []

for at in atom_type:
    atom_num.append(count_TRUE(atoms==at))
    for i in range(len(atoms)):
        if atoms[i] == at:
            pos_sort.append(pos_frac[i,:])

pos_sort = np.array(pos_sort)

lat_str = '\n'.join([' '.join([f'{item:19.8f}' for item in line]) for line in decell])
pos_str = '\n'.join([' '.join([f'{item:19.8f}' for item in line]) for line in pos_sort])
atom_str = ' '.join([f'{item:6s}' for item in atom_type])
atom_num_str = ' '.join([f'{item:6.0f}' for item in atom_num])

out_poscar = f'''sqs structure
1.0
{lat_str}
 {atom_str}
 {atom_num_str}
Direct
{pos_str}
'''

with open(f'{file_name}-POSCAR', 'w') as obj:
    obj.write(out_poscar)
