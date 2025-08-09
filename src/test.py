from adpy import run_dock

prep_ligand = '../examples/prepared_ligands/example_prepared_ligand.pdbqt'
prep_receptor = '../examples/prepared_proteins/example_prepared_protein.pdbqt'

run_dock(ligand=prep_ligand, receptor=prep_receptor)