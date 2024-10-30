# special-quasidisorder-structure-code-SQDS

This respository saves the Python scripts about the Special Quasirandom Structure (SQS), Special Quasidisorder Structure (SQDS), and alloy bandgap fitting methods.

SQS denotes that the structure is fully random in atom arrangement. Each atomic correlation function aligns with the $(2x-1)^k$, where $x$ is the alloy composition and $k$ is the atom number in the cluster.

SQDS is the structure that has a partial disorder atom arrangement, not the fully random. This code saves in `2_SQDS/`.

Disordered supercell usually has an anion-centered polyhedron distribution. The code that calculates the distribution saves in `3_Polyhedron_Distribution/`.

We newly define the disordered alloy bandgap by the electronic density of states. The code saves in `4_Fitting_Alloy_Bandgap/`. 

## 1. Converting from ATAT structure file to POSCAR file

The other converting code from ATAT structure file to POSCAR file has a atom number limition that about 400 atoms. It is hard to sufficient the large simulation.

So I write this code for no any nolimition. 

I provide the script that named `sqs2poscar-HPLiang.py`, and the example ATAT file `lat.in`.

The reader enter
```
python sqs2poscar-HPLiang.py sqs.out
```
and can obtain a POSCAR file `sqs.out-POSCAR`.

## 2. SQDS code


## 3. Polyhedron distribution code


## 4. Fitting alloy bandgap code


