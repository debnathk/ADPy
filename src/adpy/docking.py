# Run docking
from vina import Vina
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ligand', required=True, help='Set the ligand path')
parser.add_argument('--receptor', required=True, help='Set the receptor path')
parser.add_argument('--output_dir', required=True, help='Set path to save docking output')

args = parser.parse_args()

v = Vina(sf_name='vina')

# Set the receptor and ligand from PDBQT files
ligand = args.ligand
receptor = args.receptor
output_dir = args.output_dir

# ligand = './prepared_ligands/example_prepared_ligand.pdbqt'
# receptor = './prepared_proteins/example_prepared_protein.pdbqt'

v.set_receptor(receptor)
print(f'Receptor: {receptor}')

v.set_ligand_from_file(ligand)
print(f'Ligand: {ligand}')
v.compute_vina_maps(center=[-.319, 5.27, 1.59], box_size=[80, 80, 80])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose(f'{output_dir}/example_docking.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=5)
v.write_poses(f'{output_dir}/example_docking.pdbqt', n_poses=5, overwrite=True)