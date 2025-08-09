from adpy import run_dock
import argparse

# ligand = '/Users/debnathk/Documents/ADPy/examples/prepared_ligands/example_prepared_ligand.pdbqt'
# receptor = '/Users/debnathk/Documents/ADPy/examples/prepared_proteins/example_prepared_protein.pdbqt'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', required=True, help='Add Ligand Path')
    parser.add_argument('--receptor', required=True, help='Add Receptor Path')
    args = parser.parse_args()

    run_dock(ligand=args.ligand, receptor=args.receptor)