#!/bin/bash
#SBATCH -J singularity_test
#SBATCH -o singularity_test.out
#SBATCH -e singularity_test.err
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -t 0-01:20
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4000

# Singularity command line options
singularity exec --nv /n/home03/phanakata/software.simg python3 in.cyl_npt.py	
