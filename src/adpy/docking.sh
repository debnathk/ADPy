#!/bin/bash

root='/Users/debnathk/Documents/ADPy'

mkdir -p ${root}/examples/prepared_ligands
mkdir -p ${root}/examples/prepared_proteins
mkdir -p ${root}/examples/results

# Activate the conda env
eval "$(conda shell.bash hook)"
conda activate adpy

# Prepare ligands
mk_prepare_ligand.py -i ${root}/examples/ligands/example_ligand.sdf -o ${root}/examples/prepared_ligands/example_prepared_ligand.pdbqt

# Prepare targets
mk_prepare_receptor.py -i ${root}/examples/proteins/example_protein.pdb -o ${root}/examples/prepared_proteins/example_prepared_protein -p -v \
--box_size 80 80 80 --box_center -.319 5.27 1.59

# Run docking
python ${root}/src/adpy/docking.py --ligand ${root}/examples/prepared_ligands/example_prepared_ligand.pdbqt --receptor ${root}/examples/prepared_proteins/example_prepared_protein.pdbqt --output_dir ${root}/examples/results > ${root}/examples/results/example_docking_output.log 2> ${root}/examples/results/example_docking_error.log 

# Save affinities to file
python ${root}/src/adpy/utils.py --output_dir ${root}/examples/results 