# Run docking
from vina import Vina
from .utils import trimName
import argparse
import pandas as pd
import os
import csv

def dock_prep(ligand=None, receptor=None, output_dir='/Users/debnathk/Documents/ADPy/examples/results'):
    pass

def run_dock_ss(ligand=None, receptor=None, output_dir=None):
    '''
    Run docking: Single Ligand - Single Receptor

    Args:
        ligand = Prepared ligand (.pdbqt)
        receptor = Prepared receptor (.pdbqt)
        output_dir = Path to save the docking output

    '''

    v = Vina(sf_name='vina')

    # Set the receptor and ligand from PDBQT files
    # ligand = './prepared_ligands/example_prepared_ligand.pdbqt'
    # receptor = './prepared_proteins/example_prepared_protein.pdbqt'

    # Set receptor
    v.set_receptor(receptor)
    print(f'Receptor: {receptor}')

    # Set ligand
    v.set_ligand_from_file(ligand)
    print(f'Ligand: {ligand}')
    v.compute_vina_maps(center=[-.319, 5.27, 1.59], box_size=[80, 80, 80])

    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(f'{output_dir}/{trimName(ligand)}_{trimName(receptor)}.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=5)
    v.write_poses(f'{output_dir}/{trimName(ligand)}_{trimName(receptor)}.pdbqt', n_poses=5, overwrite=True)

    print(f"Docking Succesfull: Output saved in {output_dir}")

    # Save output CSV file
    output_csv = f"{trimName(ligand)}_{trimName(receptor)}_docking_results.csv"
    output_path = os.path.join(output_dir, output_csv)

    # Collect results
    results = []

    # for filename in os.listdir(results_folder):
    # if filename.endswith("output.log"):
    #     filepath = os.path.join(results_folder, filename)
    #     print(filepath)
    # ligand = receptor = affinity = None

    pdbqt_files = [f for f in os.listdir(output_dir) if f.endswith(".pdbqt")]
    pdbqt_file = os.path.join(output_dir, pdbqt_files[0])
    # print(f'Reading {pdbqt_file}')

    with open(pdbqt_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                affinity = line.split()[3]
                # print(affinity)
                break

        # Preprocess the long names of ligand, and proteins - only keep the entity name
        ligand = trimName(ligand)
        receptor = trimName(receptor)


        print(ligand, receptor, affinity)

        # if ligand and receptor and affinity:
        results.append([ligand, receptor, affinity])

        # df = pd.DataFrame({'Ligand': ligand, 'Receptor': receptor, 'Binding_Affinity(kcal/mol)': affinity}, index=[0])

    # Write to CSV
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ligand", "receptor", "binding_affinity(kcal/mol)"])
        writer.writerows(results)

    print(f"Results saved to {output_path}")

def run_dock_ms(ligand_dir=None, receptor=None, output_dir=None):
    '''
    Run docking: Single Ligand - Single Receptor

    Args:
        ligand = Prepared ligand (.pdbqt)
        receptor = Prepared receptor (.pdbqt)
        output_dir = Path to save the docking output

    '''

    v = Vina(sf_name='vina')

    # Set the receptor and ligand from PDBQT files
    # ligand = './prepared_ligands/example_prepared_ligand.pdbqt'
    # receptor = './prepared_proteins/example_prepared_protein.pdbqt'

    # Set receptor
    v.set_receptor(receptor)
    receptor = trimName(receptor)
    print(f'Receptor: {receptor}')

    # Collect results
    results_ligand = []
    results_receptor = []
    results_affinity = []

    ligands = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]
    for ligand_i in ligands:

        # Set ligand
        ligand = os.path.join(ligand_dir, ligand_i)
        v.set_ligand_from_file(ligand)
        ligand = trimName(ligand)
        print(f'Ligand: {ligand}')
        v.compute_vina_maps(center=[-.319, 5.27, 1.59], box_size=[80, 80, 80])

        # Score the current pose
        energy = v.score()
        print('Score before minimization: %.3f (kcal/mol)' % energy[0])

        # Minimized locally the current pose
        energy_minimized = v.optimize()
        print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
        # v.write_pose(f'{output_dir}/{ligand}_{receptor}.pdbqt', overwrite=True)

        # Dock the ligand
        # v.dock(exhaustiveness=32, n_poses=5)
        # v.write_poses(f'{output_dir}/{ligand}_{receptor}.pdbqt', n_poses=5, overwrite=True)

        print(f"Docking Succesfull: Output saved in {output_dir}")

        # Save ligand, receptor for later
        results_ligand.append(ligand)
        results_receptor.append(receptor)

    # Save output CSV file
    output_csv = f"{receptor}_docking_results.csv"
    output_path = os.path.join(output_dir, output_csv)

    pdbqt_files = [f for f in os.listdir(output_dir) if f.endswith(".pdbqt")]
    for pdbqt_file in pdbqt_files:
        file = os.path.join(output_dir, pdbqt_file)

        with open(file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    affinity = line.split()[3]
                    # print(affinity)
                    break


        # print(ligand, receptor, affinity)
        # if affinity:
        results_affinity.append(affinity)

    df = pd.DataFrame({'ligand': results_ligand, 'receptor': results_receptor, 'binding_affinity(kcal/mol)': results_affinity})
    print(df)
    df.to_csv(output_path, index=False)

    # Write to CSV
    # with open(output_path, "w", newline="") as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(["ligand", "receptor", "binding_affinity(kcal/mol)"])
    #     writer.writerows(results)

    print(f"Results saved to {output_path}")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--ligand', required=True, help='Set the ligand path')
#     parser.add_argument('--receptor', required=True, help='Set the receptor path')
#     parser.add_argument('--output_dir', required=True, help='Set path to save docking output')
#     args = parser.parse_args()

#     run_dock(ligand=args.ligand, receptor=args.receptor, output_dir=args.output_dir)