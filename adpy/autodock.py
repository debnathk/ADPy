import os
import csv
from typing import List, Tuple, Optional, Union
import subprocess
import requests
import sys

from vina import Vina
from .utils import trimName, extractBindingAffinity
import polars as pl


class AutoDock:
    def __init__(self, sf_name: str = "vina") -> None:
        """
        Initialize AutoDock with Vina scoring function.

        Args:
            sf_name: Scoring function name for Vina
        """
        self.v = Vina(sf_name=sf_name)
        self.default_center = (-0.319, 5.27, 1.59)
        self.default_box_size = (80, 80, 80)

    def singleLigandSingleReceptor(
        self,
        ligand: str,
        receptor: str,
        output_dir: str,
        AlphaFold: bool = True,
        center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[int, int, int]] = None,
        exhaustiveness: int = 32,
        n_poses: int = 5,
        save_csv: bool = True,
    ) -> dict:
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

        # use default values if AlphaFold protein
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

    def multiLigandSingleReceptor(
        self,
        ligand_dir: str,
        receptor: str,
        output_dir: str,
        AlphaFold: bool = True,
        center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[int, int, int]] = None,
        exhaustiveness: int = 32,
        n_poses: int = 5,
        save_csv: bool = True,
    ) -> None:
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
        ligands = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]
        df_docking_results_all = []

        try:
            for ligand in ligands:
                ligand = os.path.join(ligand_dir, ligand)
                # Setup docking
                docking_results = self._setup_and_dock(
                    ligand,
                    receptor,
                    center,
                    box_size,
                    exhaustiveness,
                    n_poses,
                    output_dir,
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

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def singleLigandMultiReceptor(
        self,
        ligand: str,
        receptor_dir: str,
        output_dir: str,
        AlphaFold: bool = True,
        center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[int, int, int]] = None,
        exhaustiveness: int = 32,
        n_poses: int = 5,
        save_csv: bool = True,
    ) -> None:
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
        receptors = [f for f in os.listdir(receptor_dir) if f.endswith(".pdbqt")]
        df_docking_results_all = []

        try:
            for receptor in receptors:
                receptor = os.path.join(receptor_dir, receptor)
                # Setup docking
                docking_results = self._setup_and_dock(
                    ligand,
                    receptor,
                    center,
                    box_size,
                    exhaustiveness,
                    n_poses,
                    output_dir,
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

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def multiLigandMultiReceptor(
        self,
        ligand_dir: str,
        receptor_dir: str,
        output_dir: str,
        AlphaFold: bool = True,
        center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[int, int, int]] = None,
        exhaustiveness: int = 32,
        n_poses: int = 5,
        save_csv: bool = True,
    ) -> None:
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
        receptors = [f for f in os.listdir(receptor_dir) if f.endswith(".pdbqt")]
        ligands = [f for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]

        df_docking_results_all = []

        try:
            for receptor in receptors:
                receptor = os.path.join(receptor_dir, receptor)
                for ligand in ligands:
                    ligand = os.path.join(ligand_dir, ligand)

                    # Setup docking
                    docking_results = self._setup_and_dock(
                        ligand,
                        receptor,
                        center,
                        box_size,
                        exhaustiveness,
                        n_poses,
                        output_dir,
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

        except Exception as e:
            print(f"Docking failed: {str(e)}")
            raise

    def _validate_input_files(self, ligand: str, receptor: str) -> None:
        """Validate that input files exist and have correct extensions.
        Run docking: Single Ligand - Single Receptor

        Args:
        ligand (str): Path to a single ligand file in .pdbqt format.
        receptor (str): Path to a single receptor file in .pdbqt format.

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        for file_path, file_type in [(ligand, "ligand"), (receptor, "receptor")]:
            if not os.path.exists(file_path):
                raise FileNotFoundError(
                    f"{file_type.capitalize()} file not found: {file_path}"
                )
            if not file_path.endswith("pdbqt"):
                raise ValueError(
                    f"{file_type.capitalize()} file must be in .pdbqt format: {file_path}"
                )

    def _validate_input_dir(self, dir: str) -> None:
        """Validate that input dir exists
        Run docking: Single Ligand - Single Receptor
        
        Args:
        input_dir (str): Path to the input directory.

        Raises:
            FileNotFoundError: If ligand or receptor file doesn't exist
            Exception: If docking process fails
        """
        if not os.path.exists(dir):
            raise FileNotFoundError(f"{dir} not found.")

    def _setup_and_dock(
        self,
        ligand: str,
        receptor: str,
        center: List[float],
        box_size: List[float],
        exhaustiveness: int,
        n_poses: int,
        output_dir: str,
    ) -> dict:
        """
    Sets up and performs molecular docking using AutoDock Vina.

    This method configures the docking parameters, scores and minimizes the 
    initial pose, performs docking, saves the resulting poses, and extracts 
    the binding affinity of the best pose.

    Args:
        ligand (str): Path to the ligand file in PDBQT format.
        receptor (str): Path to the receptor file in PDBQT format.
        center (List[float]): Coordinates [x, y, z] for the center of the docking box.
        box_size (List[float]): Dimensions [x, y, z] of the docking box.
        exhaustiveness (int): Exhaustiveness of the global search. Higher values increase accuracy and time.
        n_poses (int): Number of binding poses to generate.
        output_dir (str): Directory where docking results will be saved.

    Returns:
        dict: A dictionary containing:
            - 'ligand' (str): Base name of the ligand file.
            - 'receptor' (str): Base name of the receptor file.
            - 'binding_affinity' (float or list): Binding affinity score(s) from docking results (in kcal/mol).
    """

        # Set receptor
        self.v.set_receptor(receptor)
        print(f"Receptor: {receptor}")

        # Set ligand
        self.v.set_ligand_from_file(ligand)
        print(f"Ligand: {ligand}")

        # Configure binding site
        self.v.compute_vina_maps(center=center, box_size=box_size)

        # Score the current pose
        energy = self.v.score()
        print("Score before minimization: %.3f (kcal/mol)" % energy[0])

        # Minimized locally the current pose
        energy_minimized = self.v.optimize()
        print("Score after minimization : %.3f (kcal/mol)" % energy_minimized[0])

        # Generate output filename
        ligand_name = trimName(ligand)
        receptor_name = trimName(receptor)
        output_file = f"{output_dir}/{ligand_name}_{receptor_name}.pdbqt"

        # Save minimized pose
        self.v.write_pose(output_file, overwrite=True)

        # Dock the ligand
        self.v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        self.v.write_poses(output_file, n_poses=n_poses, overwrite=True)

        # Extract binding affinity from output
        binding_affinity = extractBindingAffinity(output_file)

        return {
            "ligand": ligand_name,
            "receptor": receptor_name,
            "binding_affinity": binding_affinity,
        }


if __name__ == "__main__":
    docker = AutoDock()
    results = docker.singleLigandSingleReceptor(
        ligand="./ligands/compound1.pdbqt",
        receptor="./receptors/protein1.pdbqt",
        output_dir="./outputs",
    )

    print(f"Binding affinity: {results['binding_affinity']} kcal/mol")


class DockPrep:
    def __init__(
        self, ligand_tool="mk_prepare_ligand.py", receptor_tool="mk_prepare_receptor.py"
    ):
        self.ligand_tool = ligand_tool
        self.receptor_tool = receptor_tool
        self.default_box_center = (-0.319, 5.27, 1.59)
        self.default_box_size = (80, 80, 80)

    def prepare_ligand(self, query: str, target: str) -> None:
        """
    Prepares a single ligand file for docking by converting or formatting it using an external tool.

    This method uses an external command-line tool (specified by `self.ligand_tool`) to process 
    the input ligand file and generate a prepared ligand file suitable for docking (e.g., in PDBQT format).

    Args:
        query (str): Path to the input ligand file (e.g., MOL2, SDF).
        target (str): Path where the prepared ligand file will be saved.

    Raises:
        Prints an error message if the external tool fails during execution.
    """
        try:
            os.makedirs(os.path.dirname(target), exist_ok=True)
            subprocess.run([self.ligand_tool, "-i", query, "-o", target], check=True)
            print(f"Ligand prepared: {target}")
        except subprocess.CalledProcessError as e:
            print(f"Error preparing ligand {query}: {e}")

    def prepare_receptor(
        self,
        query: str,
        target_prefix: str,
        AlphaFold: bool = True,
        box_size: Optional[Tuple[int, int, int]] = None,
        box_center: Optional[Tuple[float, float, float]] = None,
    ) -> None:
        """
    Prepares a single receptor file for docking by processing the input structure 
    and defining the docking box using an external tool.

    If the receptor is from AlphaFold, default box size and center (defined in 
    `self.default_box_size` and `self.default_box_center`) are used. Otherwise, 
    both `box_size` and `box_center` must be provided.

    Args:
        query (str): Path to the input receptor file (e.g., in PDB format).
        target_prefix (str): Output path prefix for the prepared receptor files.
        AlphaFold (bool, optional): Whether the receptor is from AlphaFold. If True, 
            uses default box parameters. Defaults to True.
        box_size (Optional[Tuple[int, int, int]]): Size of the docking box (x, y, z).
            Required if AlphaFold is False.
        box_center (Optional[Tuple[float, float, float]]): Center of the docking box (x, y, z).
            Required if AlphaFold is False.

    Raises:
        ValueError: If `AlphaFold` is False and either `box_size` or `box_center` is not provided.
        Prints an error message if the external tool fails during execution.
    """

        # use default values if AlphaFold protein
        if AlphaFold:
            box_center = self.default_box_center
            box_size = self.default_box_size
        else:
            if box_center is None or box_size is None:
                raise ValueError(
                    "box_size and box_center must be provided if AlphaFold=False"
                )

        try:
            os.makedirs(os.path.dirname(target_prefix), exist_ok=True)
            subprocess.run(
                [
                    self.receptor_tool,
                    "-i",
                    query,
                    "-o",
                    target_prefix,
                    "-p",
                    "-v",
                    "--box_size",
                    *map(str, box_size),
                    "--box_center",
                    *map(str, box_center),
                ],
                check=True,
            )
            print(f"Receptor prepared: {target_prefix}")
        except subprocess.CalledProcessError as e:
            print(f"Error preparing receptor {query}: {e}")

    def prepare_ligands_batch(self, ligands: List[Tuple[str, str]]) -> None:
        """
    Prepares multiple ligand files for docking in a batch process.

    This method iterates over a list of ligand input/output file pairs and prepares each 
    ligand using the `prepare_ligand` method.

    Args:
        ligands (List[Tuple[str, str]]): A list of tuples, where each tuple contains:
            - input_file (str): Path to the input ligand file.
            - output_file (str): Path where the prepared ligand file will be saved.

    Returns:
        None
    """
        for query, target in ligands:
            self.prepare_ligand(query, target)

    def prepare_receptors_batch(
        self,
        receptors: List[
            Tuple[
                str,
                str,
                bool,
                Optional[Tuple[int, int, int]],
                Optional[Tuple[float, float, float]],
            ]
        ],
    ) -> None:
        """
    Prepares multiple receptor files for docking in a batch process.

    This method iterates over a list of receptor preparation parameters and prepares each 
    receptor using the `prepare_receptor` method. If the receptor is from AlphaFold, default 
    box size and center are used.

    Args:
        receptors (List[Tuple[str, str, bool, Optional[Tuple[int, int, int]], Optional[Tuple[float, float, float]]]]): 
            A list of tuples where each tuple contains:
            - input_file (str): Path to the input receptor file.
            - output_prefix (str): Output path prefix for the prepared receptor files.
            - AlphaFold (bool): Whether the receptor is from AlphaFold.
            - box_size (Optional[Tuple[int, int, int]]): Size of the docking box (x, y, z), required if AlphaFold is False.
            - box_center (Optional[Tuple[float, float, float]]): Center of the docking box (x, y, z), required if AlphaFold is False.

    Returns:
        None
    """
        for query, target_prefix, AlphaFold, box_size, box_center in receptors:
            if AlphaFold:
                self.prepare_receptor(query, target_prefix, AlphaFold, None, None)
            else:
                self.prepare_receptor(
                    query, target_prefix, AlphaFold, box_size, box_center
                )


class getAlphaFold:
    def __init__(self, gene: Optional[str] = None, uniprot: Optional[str] = None) -> None:
        self.gene = gene
        self.uniprot = uniprot
    
    def _from_gene_single(self, gene: str):
        params = {
            "query": f"{gene} AND reviewed:true AND gene_exact:{gene}",
            "fields": ["accession", "protein_name"],
            #   "sort": "accession desc",
            "size": "50",
        }
        headers = {"accept": "application/json"}
        base_url = "https://rest.uniprot.org/uniprotkb/search?query=(organism_id:9606)"
        
        response = requests.get(base_url, headers=headers, params=params)
        if not response.ok:
            response.raise_for_status()
            sys.exit()
        try:
            data = response.json()
            accs = data["results"][0]["primaryAccession"]
            
        except Exception as e:
            print(f"Error occured for {gene}: {str(e)}")
            return None
        
        print(f"Gene: {gene} --> UniProt ID: {accs}")
        self._from_uniprot_single(accs)
    
    def _from_uniprot_single(self, uniprot: str) -> None:
        af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
        response = requests.get(af_url)
        # if not response.ok:
        #     response.raise_for_status()
        #     sys.exit()
        try:
            data = response.json()
            af_id = data[0]["entryId"]
            
            print(f"UniProt ID: {uniprot} --> AlphaFold ID: {af_id}")
            
            pdb_url = data[0]["pdbUrl"]
            
            pdb_path = './examples/receptors/'
            os.makedirs(pdb_path, exist_ok=True)
            output_path = os.path.join(pdb_path, f"{af_id}.pdb")
            
            pdb_response = requests.get(pdb_url)
            if pdb_response.ok:
                with open(output_path, "wb") as f:
                    f.write(pdb_response.content)
                print(f"Downloaded PDB file as {output_path}")
            else:
                print(f"Error while downloading PDB for {af_id}. HTTPS error: {pdb_response.status_code}")
        except Exception as e:
            print(f"Error occured while fetching data for {uniprot}: {str(e)}")

    def _from_gene_batch(self, gene_list: List[str]) -> None:
        for gene in gene_list:
            self._from_gene_single(gene)

    def _from_uniprot_batch(self, uniprot_list: List[str]) -> None:
        for accs in uniprot_list:
            self._from_uniprot_single(accs)
