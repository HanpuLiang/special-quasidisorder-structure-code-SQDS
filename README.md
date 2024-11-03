# special-quasidisorder-structure-code-SQDS

This respository saves the Python scripts about the Special Quasirandom Structure (SQS), Special Quasidisorder Structure (SQDS), and alloy bandgap fitting methods.

SQS denotes that the structure is fully random in atom arrangement. Each atomic correlation function aligns with the $(2x-1)^k$, where $x$ is the alloy composition and $k$ is the atom number in the cluster.

SQDS is the structure that has a partial disorder atom arrangement, not only the fully random. This code saves in `2_SQDS/`.

Disordered supercell usually has an anion-centered polyhedron distribution. The code that calculates the distribution saves in `3_Polyhedron_Distribution/`.

We newly define the disordered alloy bandgap by the electronic density of states. The code saves in `4_Fitting_Alloy_Bandgap/`. 

**These codes are applied into my papers: [JACS 146, 16222 (2024)](https://pubs.acs.org/doi/10.1021/jacs.4c04201), so I hope that if the reader uses my code, please cite these paper, especially the SQDS code. Thanks!**

bib content:
```

```

## 1. Converting from ATAT structure file to POSCAR file

The other converting code from ATAT structure file to POSCAR file has a atom number limition that about 400 atoms. It is hard to sufficient the large simulation.

So I write this code for no any nolimition. 

I provide the script that named `sqs2poscar-HPLiang.py`, and the example ATAT file `sqs.out`.

The reader uses this code by entering
```linux
python sqs2poscar-HPLiang.py sqs.out
```
and can obtain a POSCAR file `sqs.out-POSCAR`.

This script is based on `Python_3.x` and requires `numpy` package.

## 2. SQDS code

This code is based on the `corrdump` command in `ATAT` code to obtain the atomic correlation function, so when the reader runs this code, the `ATAT` must be installed in your computer/cluster.

The original `mcsqs` command in `ATAT` is hard to generate a partial disordered structure and I don't know why it can not work effectively. So I write this SQDS to construct the partial disordered structure that satisfy the targeted atomic correlation function.

Furthermore, `mcsqs` also don't consider the matching of long-range order, which has another definition: 
```
In the ordered structure, A and B seperate into two sublattices. When the structure becomes more disorder, the A and B continus to exchange between both sublattices, while A and B are uniformly distribute in every sublattice. At the fully disordered state, A and B are half half seperate into two sublattices and in each sublattice, the atoms distributions are random. 
```
`mcsqs` only search the structures by the targeted random correlation function, while ignore the exchange between two sublattice from the ordered structure. Thus, although the short-range order, i.e., atomic correlation, closes to 0 of this searched SQS, it may is not satisfy the long-range order closes to 0. So, this SQDS code search disordered structure not only matches the correlation function, but also satisfies the atom exchange ratio of long-range order.

In SQDS method, the correlation function at targeted long-range order $\eta$ is defined as
$$
\Pi(x, \eta) = \Pi(x,0)+ \eta^2[\Pi(x,1)-\Pi(x,0)],
$$
where $\Pi(x,0)=(2x-1)^k$ and $\Pi(x,1)$ is the correlation function of the targeted ordered structure, such as the chalcopyrite of ZnSnP2 compound.

We set some parameters in the front of code, such as
```python
target_lro = 0.0    # the targeted Long-range order $\eta$
target_sro = 0.0    # the targeted short-range order, i.e., the atomic correlation function of each cluster
N_iter = 100000     # iteration number
cutoff_score = 0.001    # the score cutoff during the search process
expand_matrix = np.array([[5, 0, 0], [0, 5, 0], [0, 0, 2]]) # the expand coefficient of targeted ordered structure
exchange_atoms = ['Zn', 'Sn']   # the exchange atoms
fixed_atoms = ['P'] # the unexchange atom
order_file_name = 'ZnSnP2-symmetry-cell.vasp'   # the filename of ordered structure
SRO_critic = 2  # 0: mean, 1: max, 2: exp(-(d-d0))
```

In the search process, the score is defined as $score=0.5*|\eta_{target} - \eta_{current}| + 0.5*|\Pi_{target} - \Pi_{current}|$. The curoff of score presents that the code only outputs the searched structure below this cutoff value. Generally, for the different lattice, composition, and order degree, we need to choose a suitable small cutoff, otherwise the code will generates a hugh number of structures or no structure generated. The reader should test this cutoff for a targeted case.

This code calculates the long-range order is calculated by the initial configuration of inputed ordered structure. I summurizes the concentration that A atom still stay in the A sublattice, and $\eta=2n_A-1$. The short-range order, i.e., atomic correlation function, is calculated by the `corrdump` command in `ATAT` code. So every generated structure will convert to the ATAT structure files and output the correlation function with the corresponding input file `clusters.out` and `rndstr.in`.

This code supports multiprocess that submitting to the computer cluster. After runing the code, many files will generate and the searched structure sill save into the `save-best-data/`. The score more smaller, the structure more closer to the target disorder degree. 

## 3. Polyhedron distribution code

In disordered ABC2 compound, while A and B are the exchangeable cations, the appearence frequency of anion-centered polyhedra $A_nB_{N-n}$ is $C^n_N$. In zinc-blende lattice, the nearest-neighbor polyhedron is tetrahedron, and the ordered configuration usually is constructed by $A_2B_2$ tetrahedron. As the structure closes to the disordered configuration, all motifs of the tetrahedra are $A_0B_4$, $A_1B_3$, $A_2B_2$, $A_3B_1$, abd $A_4B_0$ with the ratio 1:4:6:4:1. The second-nearest-neghbor polyhedron has 12 vertices.

This code takes zinc-blende ZnSnP2 alloys as the example, summarzing all nearest-neighbor and second-nearest-neghbor anion-centered polyhedra in the input structure. 

The reader can motify the `mats` list and `cation`, `anion` variables in the code and run the code by
```
python polyhedron_distribution.py
```

This code is dependented on the `ase` package.

## 4. Fitting alloy bandgap code

This method is proposed in the paper `xxxx` (ready to submit). In real physics, the bandgap of disordered alloy usually has a measureable value by optical measurement. However, the bandgap usually is hard to converge with the simulated supercell increasing. This is because the size-dependent localized states in the band gap. So we solve this problem by a fitting method based on the electronic density of states.

We assume that the disordered alloys have the nonparabolic band near the band edge with the relation
$\gamma(E)=E+c_1E^2=\hbar k^2/2/m_d$. Accourding to the relation between density of states and $\gamma(E)$, $N(E)\sim \gamma^{0.5}(d\gamma/dE)$, we can obtain the nonparabolic density of states is $N^2(E)=E+5c_1E^2+8c_2^2E^3+4c_1^3E^4$. Thereby, we can fit the band edge states by this formula and redefine the VBM and CBM by the extrapolation value of fitting line at $N^2(E)=0$.

As the case of the ZnSnP2 in 64-atom and 512-atom supercells, the corrected bandgaps both converge to the 1.0 eV, while the occupation bandgap is 0.6 and 0.1 eV, respectively.

I provide two versions of this code, Python and Matlab, to get the results in the paper. The reader can directly run the code to get the results.

Mention that, for a new DOS in a new compounds, it is important to choose a suitable fitting range to obtain a reasonable fitting results. The choose of the disordered supercell also decide the results. For example, if a supercell only satisfies the short-range order but ignore the long-range order, the DOS may be not very useful. Before we calculate the bandgap, we must to ensure that the supercell is correct.
