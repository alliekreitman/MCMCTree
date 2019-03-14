#!/bin/bash

#SBATCH --job-name=MCMCTree_p1
#SBATCH --error=Tree.txt
#SBATCH --time=24:00:00
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --partition=sandyb
#SBATCH --mem=20g

module load openmpi

mpiexec ./out
