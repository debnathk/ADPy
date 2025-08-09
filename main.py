from adpy import run_dock, pdbqt2csv
import argparse
import os

# ligand = '/Users/debnathk/Documents/ADPy/examples/prepared_ligands/example_prepared_ligand.pdbqt'
# receptor = '/Users/debnathk/Documents/ADPy/examples/prepared_proteins/example_prepared_protein.pdbqt'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', required=True, help='Set ligand path')
    parser.add_argument('--receptor', required=True, help='Set receptor path')
    parser.add_argument('--output-dir', required=True, help='Set path to save docking output')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    run_dock(ligand=args.ligand, receptor=args.receptor, output_dir=args.output_dir)
    # pdbqt2csv(ligand=args.ligand, receptor=args.receptor, output_dir=args.output_dir)
    