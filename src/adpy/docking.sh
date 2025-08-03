
# Prepare ligands
mk_prepare_ligand.py -i ./ligands/Structure2D_COMPOUND_CID_24978538.sdf -o ./prepared_ligands/ligand_24978538.pdbqt

# Prepare targets
mk_prepare_receptor.py -i ./proteins/AF-A5LHX3-F1.pdb -o ./prepared_proteins/AF-A5LHX3-F1_receptor -p -v \
--box_size 80 80 80 --box_center -.3190 5.2703 1.59

# Run docking
from vina import Vina


v = Vina(sf_name='vina')

# Set the receptor and ligand from PDBQT files
ligand = './prepared_ligands/ligand_24978538.pdbqt'
receptor = './prepared_proteins/AF-A5LHX3-F1_receptor.pdbqt'

v.set_receptor(receptor)
print(f'Receptor: {receptor}')

v.set_ligand_from_file(ligand)
print(f'Ligand: {ligand}')
v.compute_vina_maps(center=[-.319, 5.27, 1.59], box_size=[67, 50, 80])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=5)
v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)