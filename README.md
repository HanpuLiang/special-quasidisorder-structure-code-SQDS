# special-quasidisorder-structure-code-SQDS

This repository contains Python scripts related to the Special Quasirandom Structure (SQS), Special Quasidisorder Structure (SQDS), and alloy bandgap fitting methods.

- SQS refers to a structure with a fully random atomic arrangement, where the atomic correlation function follows $(2x-1)^k$, where $x$ is the alloy composition and $k$ is the number of atoms in the cluster.
- SQDS represents a structure with a partially disordered atomic arrangement, as opposed to fully random. The related code is located in the `2_SQDS/` directory.
- A disordered supercell typically exhibits an anion-centered polyhedron distribution. The script for calculating the distribution is located in the `3_Polyhedron_Distribution/` directory.
- A new method for defining the disordered alloy bandgap based on the electronic density of states is implemented. The related code can be found in `4_Fitting_Alloy_Bandgap/` directory.

**These scripts are used in my paper: [JACS 146, 16222 (2024)](https://pubs.acs.org/doi/10.1021/jacs.4c04201). Another paper is now being submitted to the PRL.**



## 1. Converting from ATAT structure file to POSCAR file

The existing code for converting ATAT structure files to POSCAR has a limitation of about 400 atoms, making it difficult to handle large simulations. I have written a new script that removes this limitation.

You can use the script `sqs2poscar-HPLiang.py` with the example ATAT file `sqs.out`. To run it, enter the following command:
```bash
python sqs2poscar-HPLiang.py sqs.out
```
This will generate a POSCAR file named `sqs.out-POSCAR`.

The script requires `Python 3.x` and the `numpy` package.

## 2. SQDS code

This code is based on the `corrdump` command from the `ATAT` code to obtain the atomic correlation function. To run this, you must have `ATAT` installed on your computer or cluster.

The original `mcsqs` command in `ATAT` is not capable of generating a partially disordered structure effectively, and I am unsure why it does not work as expected. Therefore, I wrote this SQDS code to construct a partially disordered structure that satisfies the targeted atomic correlation function.

Additionally, the mcsqs command does not consider long-range order, which has an alternative definition:

> In an ordered structure, atoms A and B form two sublattices. As the structure becomes more disordered, atoms A and B exchange between the sublattices, eventually resulting in a uniform distribution. In a fully disordered state, A and B are evenly distributed across both sublattices, and the atoms within each sublattice are randomly distributed.

`mcsqs` searches for structures based solely on the short-range order (atomic correlation function) and ignores the long-range order (the exchange between sublattices). As a result, although the short-range order approaches zero in the generated SQS structure, the long-range order may still deviate from zero.

The SQDS code, therefore, searches for disordered structures that both match the correlation function and satisfy the exchange ratio of the long-range order.

### Long-range order definition

In the SQDS method, the correlation function at the targeted long-range order $\eta$ is defined as: $\Pi(x, \eta) = \Pi(x,0)+ \eta^2[\Pi(x,1)-\Pi(x,0)]$, where $\Pi(x,0)=(2x-1)^k$ and $\Pi(x,1)$ is the correlation function of the targeted ordered structure, such as the chalcopyrite structure of ZnSnP2.

### Parameters in the code

The following parameters are defined at the beginning of the code:
```python
target_lro = 0.0    # targeted Long-range order $\eta$
target_sro = 0.0    # targeted short-range order, i.e., atomic correlation function of each cluster
N_iter = 100000     # iteration number
cutoff_score = 0.001    # score cutoff during search process
expand_matrix = np.array([[5, 0, 0], [0, 5, 0], [0, 0, 2]]) # expand coefficient of targeted ordered structure
exchange_atoms = ['Zn', 'Sn']   # changeable atoms
fixed_atoms = ['P'] # fixed atom
order_file_name = 'ZnSnP2-symmetry-cell.vasp'   # filename of ordered structure
SRO_critic = 2  # 0: mean, 1: max, 2: exp(-(d-d0))
```

During the search, the score is defined as: $score=0.5*|\eta_{target} - \eta_{current}| + 0.5*|\Pi_{target} - \Pi_{current}|$. The cutoff value filters the structures based on this score. A smaller cutoff produces structures closer to the target disorder degree. For different lattices, compositions, and order degrees, the appropriate cutoff should be tested.

The long-range order is calculated from the initial configuration of the input ordered structure, where the concentration of A atoms remaining in the A sublattice is summarized, and $\eta=2n_A-1$. The short-range order is calculated using the `corrdump` command from `ATAT`. Thus, each generated structure is converted to ATAT structure files, with corresponding correlation functions and input files `clusters.out` and `rndstr.in`.

This code supports multiprocessing for cluster use. After running, the generated structures will be saved in the `save-best-data/` directory, and the structures with the smallest score will be output.

## 3. Polyhedron distribution code

In disordered ABC2 compound, while A and B are the exchangeable cations, the appearence frequency of anion-centered polyhedra $A_nB_{N-n}$ follows the combinatorial coefficients $C^n_N$. In zinc-blende lattice, the nearest-neighbor polyhedron is a tetrahedron, and the ordered configuration is usually constructed with $A_2B_2$ tetrahedron. As the structure becomes more disordered, the tetrahedral motifs vary, including $A_0B_4$, $A_1B_3$, $A_2B_2$, $A_3B_1$, abd $A_4B_0$ with a ratio 1:4:6:4:1. The second-nearest-neghbor polyhedron has 12 vertices.

This code takes zinc-blende ZnSnP2 alloys as an example and calculates the nearest-neighbor and second-nearest-neighbor anion-centered polyhedra in the input structure.

To run the code, modify the `mats` list and the `cation` and `anion` variables in the script, and then execute:
```bash
python polyhedron_distribution.py
```

This script depends on the `ase` package.

## 4. Fitting alloy bandgap code

This method, proposed in an upcoming paper, addresses the challenge of calculating the bandgap of disordered alloys. Due to size-dependent localized states near the band edge, the bandgap does not always converge with increasing supercell size. To solve this, we propose a fitting method based on the electronic density of states (DOS).

We assume the disordered alloys exhibit a nonparabolic band near the band edge, with the relation: $\gamma(E)=E+c_1E^2=\hbar k^2/2/m_d$. From the relation between DOS and $\gamma(E)$, $N(E)\sim \gamma^{0.5}(d\gamma/dE)$, we derive a fitting formula for the nonparabolic DOS: $N^2(E)=E+5c_1E^2+8c_2^2E^3+4c_1^3E^4$. This formula allows us to fit the band-edge states and extrapolate the VBM and CBM by finding the value of $E$ where $N^2(E)=0$.

For ZnSnP2 in 64-atom and 512-atom supercells, the corrected bandgaps both converge to 1.0 eV, while the uncorrected bandgaps are 0.6 and 0.1 eV, respectively.

I provide both Python and Matlab versions of this code. The reader can run the code to obtain the results presented in the paper.

Note: For new compounds, it is important to select an appropriate fitting range for reliable results. The choice of disordered supercell also influences the accuracy of the DOS. If the supercell only satisfies short-range order and ignores long-range order, the DOS may not be very useful. Therefore, ensure that the supercell is correct before calculating the bandgap.