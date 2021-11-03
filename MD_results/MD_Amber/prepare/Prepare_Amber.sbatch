#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=MD_prep

module purge
module load pdb2pqr/3.1.0
module load openbabel/intel/3.1.1
module load amber/openmpi/intel/20.06

lig='13013_2_6vzk'
pro='6vzk'

## add hydrogens for protein using pdb2pqr, pH=7.0
pdb2pqr30 --ff=AMBER --ffout=AMBER --titration-state-method=propka --with-ph=7.0 \
  ${pro}.pdb ${pro}.pqr

## add hydrogens for docked ligand using openbabel
obabel -ipdb ${lig}.pdb -opdb -O UNL.pdb --AddNonPolarH
grep -E 'HETATM |ATOM ' UNL.pdb > UNL_new.pdb

## prepare force field parameters for docked ligand
antechamber -i UNL_new.pdb -fi pdb -o UNL.mol2 -fo mol2 -c bcc -s 2
parmchk2 -i UNL.mol2 -f mol2 -o UNL.frcmod
tleap -f tleap1.in
## remove files
rm ANTECHAMBER* ATOMTYPE* sqm*

## prepare prmtop and inpcrd files for MD system (Amber)
sed -i "s/6vzk/$pro/g" tleap2.in
tleap -f tleap2.in


