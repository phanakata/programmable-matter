#!/bin/bash
#SBATCH -J cylHOOMD
#SBATCH -o HOOMD_cyl.out
#SBATCH -e HOOMD_cyl.err
#SBATCH -p shared
#SBATCH -t 0-01:30
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4000

# Singularity command line options
singularity exec /n/home03/phanakata/software.simg python3 in.cyl_npt.py
