#!/bin/bash

mkdir -p ./examples/prepared_ligands
mkdir -p ./examples/prepared_proteins
mkdir -p ./examples/results

# Prepare ligands
mk_prepare_ligand.py -i ./examples/ligands/example_ligand.sdf -o ./examples/prepared_ligands/example_prepared_ligand.pdbqt

# Prepare targets
mk_prepare_receptor.py -i ./examples/proteins/example_protein.pdb -o ./examples/prepared_proteins/example_prepared_protein.pdbqt -p -v \
--box_size 80 80 80 --box_center -.319 5.27 1.59

# Run docking
python docking.py --ligand ./examples/prepared_ligands/example_prepared_ligand.pdbqt --receptor ./examples/prepared_proteins/example_prepared_protein.pdbqt