#!/bin/bash

#SBATCH --job-name=MD
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=4:00:00
#SBATCH --mem=4GB


module purge
module load amber/openmpi/intel/20.06

mpirun -np 16 pmemd.MPI -O -i ../heat.in -o heat.out -p ../com_solvated.prmtop \
	-c ../min3.rst -r heat.rst -x heat.mdcrd -ref ../min3.rst
mpirun -np 16 pmemd.MPI -O -i ../density.in -o density.out -p ../com_solvated.prmtop \
	-c heat.rst -r density.rst -x density.mdcrd -ref heat.rst
mpirun -np 16 pmemd.MPI -O -i ../equil.in -o equil.out -p ../com_solvated.prmtop \
	-c density.rst -r equil.rst -x equil.mdcrd -ref density.rst


