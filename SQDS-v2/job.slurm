#!/bin/bash
#SBATCH -J SQDS
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -p cp6
#SBATCH -N 1
#SBATCH -n 56
cd $SLURM_SUBMIT_DIR

module load vasp/5.4.4

yhrun -n 56 python SQDS-LRO-SRO-HPLiang.py > log
