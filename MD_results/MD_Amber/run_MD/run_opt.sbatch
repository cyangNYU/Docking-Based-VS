#!/bin/bash

#SBATCH --job-name=MD
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --time=2:00:00
#SBATCH --mem=4GB


module purge
module load amber/openmpi/intel/20.06

mpirun -np 8 sander.MPI -O -i min1.in -o min1.out -p com_solvated.prmtop \
       	-c com_solvated.inpcrd -r min1.rst -ref com_solvated.inpcrd
mpirun -np 8 sander.MPI -O -i min2.in -o min2.out -p com_solvated.prmtop \
	-c min1.rst -r min2.rst -ref min1.rst
mpirun -np 8 sander.MPI -O -i min3.in -o min3.out -p com_solvated.prmtop \
	-c min2.rst -r min3.rst -ref min2.rst


