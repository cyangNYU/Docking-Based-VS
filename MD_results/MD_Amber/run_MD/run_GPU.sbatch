#!/bin/bash

#SBATCH --job-name=MD
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --tasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=4GB


module purge
module load amber/openmpi/intel/20.06

pmemd.cuda -O -i ../prod.in -o prod.out -p ../com_solvated.prmtop \
	-c equil.rst -r prod.rst -x prod.mdcrd



