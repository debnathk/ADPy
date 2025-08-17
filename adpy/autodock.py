import os
import csv
from typing import List, Tuple, Optional, Union

from vina import Vina
from .utils import trimName, extractBindingAffinity
import polars as pl

class AutoDock:
    def __init__(self, sf_name: str = 'vina') -> None:
        """
        Initialize AutoDock with Vina scoring function.

        Args:
            sf_name: Scoring function name for Vina
        """
        self.v = Vina(sf_name=sf_name)
        self.default_center = [-.319, 5.27, 1.59]
        self.default_box_size = [80, 80, 80]

    def singleLigandSingleReceptor(self, ligand: str, receptor: str, output_dir: str, AlphaFold: bool = True, center: Optional[List[float]] = None, box_size: Optional[List[float]] = None, exhaustiveness: int = 32, n_poses: int = 5, save_csv: bool = True) -> dict:
        """
        Run docking: Single Ligand - Single Receptor

        Args:
            ligand: Path to prepared ligand (.pdbqt)
            receptor: Path to prepared receptor (.pdbqt)
            output_dir: Directory to save docking output
            centre: pass
            box_size: pass
            exhaustiveness: pass
            n_poses: pass
            save_csv: pass

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        # Validate input files
        self._validate_input_files(ligand, receptor)

        # use default values if not provided
        center = self.default_center if AlphaFold else center
        box_size = self.default_box_size if AlphaFold else box_size

        # Ensure output directories exist
        os.makedirs(output_dir, exist_ok=True)

        try:
            # Setup docking
            docking_results = self._setup_and_dock(
                ligand, receptor, center, box_size, exhaustiveness, n_poses, output_dir
            )

            # Save results if requested
            if save_csv:
                # self._save_results_to_csv(docking_results, output_dir)
                df_docking_results = pl.DataFrame(docking_results)
                csv_path = f"{docking_results['ligand']}_{docking_results['receptor']}_docking_results.csv"
                output_path = os.path.join(output_dir, csv_path)
                df_docking_results.write_csv(output_path)

            print(f"Docking successful: Output saved in {output_dir}")
            return docking_results

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise 

    def multiLigandSingleReceptor(self, ligand_dir: str, receptor: str, output_dir: str, AlphaFold: bool = True, center: Optional[List[float]] = None, box_size: Optional[List[float]] = None, exhaustiveness: int = 32, n_poses: int = 5, save_csv: bool = True) -> None:
        """
        Run docking: Single Ligand - Single Receptor

        Args:
            ligand: Path to prepared ligand (.pdbqt)
            receptor: Path to prepared receptor (.pdbqt)
            output_dir: Directory to save docking output
            centre: pass
            box_size: pass
            exhaustiveness: pass
            n_poses: pass
            save_csv: pass

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        # Validate input files
        self._validate_input_dir(ligand_dir)

        # use default values if not provided
        center = self.default_center if AlphaFold else center
        box_size = self.default_box_size if AlphaFold else box_size

        # Ensure output directories exist
        os.makedirs(output_dir, exist_ok=True)

        # List all the ligands from the directory
        ligands = [f for f in os.listdir(ligand_dir) if f.endswith('.pdbqt')]
        df_docking_results_all = []

        try:
            for ligand in ligands:
                ligand = os.path.join(ligand_dir, ligand)
                # Setup docking
                docking_results = self._setup_and_dock(
                    ligand, receptor, center, box_size, exhaustiveness, n_poses, output_dir
                )
                df_docking_results = pl.DataFrame(docking_results)
                df_docking_results_all.append(df_docking_results)

            df_final = pl.concat(df_docking_results_all)
            output_csv = f"{df_final['receptor'][0]}_docking_results.csv"
            output_path = os.path.join(output_dir, output_csv)
            df_final.write_csv(output_path)

            # Check if results not saved as CSV
            if not os.path.exists(output_path):
                print(f"Failed to create {output_path}")
            print(f"Docking successful: Output saved in {output_dir}")
            # # Save results if requested
            # if save_csv:
            #     self._save_results_to_csv(df_final, output_dir)

            # return docking_results

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def singleLigandMultiReceptor(self, ligand: str, receptor_dir: str, output_dir: str, AlphaFold: bool = True, center: Optional[List[float]] = None, box_size: Optional[List[float]] = None, exhaustiveness: int = 32, n_poses: int = 5, save_csv: bool = True) -> None:
        """
        Run docking: Single Ligand - Single Receptor

        Args:
            ligand: Path to prepared ligand (.pdbqt)
            receptor: Path to prepared receptor (.pdbqt)
            output_dir: Directory to save docking output
            centre: pass
            box_size: pass
            exhaustiveness: pass
            n_poses: pass
            save_csv: pass

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        # Validate input files
        self._validate_input_dir(receptor_dir)

        # use default values if not provided
        center = self.default_center if AlphaFold else center
        box_size = self.default_box_size if AlphaFold else box_size

        # Ensure output directories exist
        os.makedirs(output_dir, exist_ok=True)

        # List all the receptors from the directory
        receptors = [f for f in os.listdir(receptor_dir) if f.endswith('.pdbqt')]
        df_docking_results_all = []

        try:
            for receptor in receptors:
                receptor = os.path.join(receptor_dir, receptor)
                # Setup docking
                docking_results = self._setup_and_dock(
                    ligand, receptor, center, box_size, exhaustiveness, n_poses, output_dir
                )
                df_docking_results = pl.DataFrame(docking_results)
                df_docking_results_all.append(df_docking_results)

            df_final = pl.concat(df_docking_results_all)
            output_csv = f"{df_final['ligand'][0]}_docking_results.csv"
            output_path = os.path.join(output_dir, output_csv)
            df_final.write_csv(output_path)

            # Check if results not saved as CSV
            if not os.path.exists(output_path):
                print(f"Failed to create {output_path}")
            print(f"Docking successful: Output saved in {output_dir}")
            # # Save results if requested
            # if save_csv:
            #     self._save_results_to_csv(df_final, output_dir)

            # return docking_results

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def multiLigandMultiReceptor(self, ligand_dir: str, receptor_dir: str, output_dir: str, AlphaFold: bool = True, center: Optional[List[float]] = None, box_size: Optional[List[float]] = None, exhaustiveness: int = 32, n_poses: int = 5, save_csv: bool = True) -> None:
        """
        Run docking: Single Ligand - Single Receptor

        Args:
            ligand: Path to prepared ligand (.pdbqt)
            receptor: Path to prepared receptor (.pdbqt)
            output_dir: Directory to save docking output
            centre: pass
            box_size: pass
            exhaustiveness: pass
            n_poses: pass
            save_csv: pass

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        # Validate ligand and receptor dirs
        self._validate_input_dir(ligand_dir)
        self._validate_input_dir(receptor_dir)

        # use default values if not provided
        center = self.default_center if AlphaFold else center
        box_size = self.default_box_size if AlphaFold else box_size

        # Ensure output directories exist
        os.makedirs(output_dir, exist_ok=True)

        # List all ligands and receptors from the directory
        receptors = [f for f in os.listdir(receptor_dir) if f.endswith('.pdbqt')]
        ligands = [f for f in os.listdir(ligand_dir) if f.endswith('.pdbqt')]

        df_docking_results_all = []

        try:
            for receptor in receptors:
                receptor = os.path.join(receptor_dir, receptor)
                for ligand in ligands:
                    ligand = os.path.join(ligand_dir, ligand)

                    # Setup docking
                    docking_results = self._setup_and_dock(
                        ligand, receptor, center, box_size, exhaustiveness, n_poses, output_dir
                    )
                    df_docking_results = pl.DataFrame(docking_results)
                    df_docking_results_all.append(df_docking_results)

            df_final = pl.concat(df_docking_results_all)
            output_csv = "docking_results.csv"
            output_path = os.path.join(output_dir, output_csv)
            df_final.write_csv(output_path)

            # Check if results not saved as CSV
            if not os.path.exists(output_path):
                print(f"Failed to create {output_path}")
            print(f"Docking successful: Output saved in {output_dir}")
            # # Save results if requested
            # if save_csv:
            #     self._save_results_to_csv(df_final, output_dir)

            # return docking_results

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def _validate_input_files(self, ligand: str, receptor: str) -> None:
        """Validate that input files exist and have correct extensions."""
        for file_path, file_type in [(ligand, "ligand"), (receptor, "receptor")]:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"{file_type.capitalize()} file not found: {file_path}")
            if not file_path.endswith('pdbqt'):
                raise ValueError(f"{file_type.capitalize()} file must be in .pdbqt format: {file_path}")
            
    def _validate_input_dir(self, dir: str) -> None:
        """Validate that input dir exists"""
        if not os.path.exists(dir):
            raise FileNotFoundError(f"{dir} not found.")
            
    def _setup_and_dock(self, ligand: str, receptor: str, center: List[float], box_size: List[float], exhaustiveness: int, n_poses: int, output_dir: str) -> dict:
        """Setup Vina parameters and perform docking."""

        # Set receptor
        self.v.set_receptor(receptor)
        print(f'Receptor: {receptor}')

        # Set ligand
        self.v.set_ligand_from_file(ligand)
        print(f'Ligand: {ligand}')

        # Configure binding site
        self.v.compute_vina_maps(center=center, box_size=box_size)

        # Score the current pose
        energy = self.v.score()
        print('Score before minimization: %.3f (kcal/mol)' % energy[0])

        # Minimized locally the current pose
        energy_minimized = self.v.optimize()
        print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])

        # Generate output filename
        ligand_name = trimName(ligand)
        receptor_name = trimName(receptor)
        output_file = f'{output_dir}/{ligand_name}_{receptor_name}.pdbqt'
        
        # Save minimized pose
        self.v.write_pose(output_file, overwrite=True)

        # Dock the ligand
        self.v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        self.v.write_poses(output_file, n_poses=n_poses, overwrite=True)

        # Extract binding affinity from output
        binding_affinity = extractBindingAffinity(output_file)

        return {
            'ligand': ligand_name,
            'receptor': receptor_name,
            'binding_affinity': binding_affinity
            # 'output_file': output_file
        }

    def _save_results_to_csv(self, results: dict, output_dir: str) -> None:
        """Save docking results to CSV file"""
        # Save output CSV file
        output_csv = f"{results['ligand']}_{results['receptor']}_docking_results.csv"
        output_path = os.path.join(output_dir, output_csv)

        # Collect results
        # df_results = pl.DataFrame(results)

        # # for filename in os.listdir(results_folder):
        # # if filename.endswith("output.log"):
        # #     filepath = os.path.join(results_folder, filename)
        # #     print(filepath)
        # # ligand = receptor = affinity = None

        # pdbqt_files = [f for f in os.listdir(output_dir) if f.endswith(".pdbqt")]
        # pdbqt_file = os.path.join(output_dir, pdbqt_files[0])
        # # print(f'Reading {pdbqt_file}')

        # with open(pdbqt_file, "r", encoding="utf-8", errors="ignore") as f:
        #     for line in f:
        #         if line.startswith("REMARK VINA RESULT:"):
        #             affinity = line.split()[3]
        #             # print(affinity)
        #             break

        #     # Preprocess the long names of ligand, and proteins - only keep the entity name
        #     ligand = trimName(ligand)
        #     receptor = trimName(receptor)


        #     print(ligand, receptor, affinity)

        #     # if ligand and receptor and affinity:
        #     results.append([ligand, receptor, affinity])

        #     # df = pd.DataFrame({'Ligand': ligand, 'Receptor': receptor, 'Binding_Affinity(kcal/mol)': affinity}, index=[0])

        # Write to CSV
        with open(output_path, "w") as csvfile:
            df_results = pl.DataFrame(results)
            df_results.write_csv(csvfile)
            # writer = csv.writer(csvfile)
            # writer.writerow(["ligand", "receptor", "binding_affinity(kcal/mol)"])
            # writer.writerow([results['ligand'], results['receptor'], results['binding_affinity']])

        print(f"Results saved to {output_path}")


if __name__ == "__main__":
    docker = AutoDock()
    results = docker.singleLigandSingleReceptor(
        ligand="./ligands/compound1.pdbqt",
        receptor="./receptors/protein1.pdbqt", 
        output_dir="./outputs",
    )

    print(f"Binding affinity: {results['binding_affinity']} kcal/mol")
