from adpy import run_dock_ss, run_dock_ms, pdbqt2csv
import argparse
import os

# ligand = '/Users/debnathk/Documents/ADPy/examples/prepared_ligands/example_prepared_ligand.pdbqt'
# receptor = '/Users/debnathk/Documents/ADPy/examples/prepared_proteins/example_prepared_protein.pdbqt'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', help='Set ligand path')
    parser.add_argument('--ligand-dir', help='Set ligand directory')
    parser.add_argument('--receptor', help='Set receptor path')
    parser.add_argument('--receptor-dir', help='Set receptor directory')
    parser.add_argument('--output-dir', required=True, help='Set path to save docking output')
    parser.add_argument('--multi-ligand', type=bool, default=False, help='If docking with multiple proteins: True, else: False')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # run_dock_ss(ligand=args.ligand, receptor=args.receptor, output_dir=args.output_dir) # one ligand - one receptor
    run_dock_ms(ligand_dir=args.ligand_dir, receptor=args.receptor, output_dir=args.output_dir) # many ligand - many receptor
    