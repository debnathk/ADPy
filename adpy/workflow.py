from typing import List, Tuple, Optional, Union
from .alphafold import AlphaFold
from .autodock import AutoDock
from .dockprep import DockPrep
from .utils import trimName, extractBindingAffinity
import polars as pl
import os

class Workflows:

    def __init__(self) -> None:
        pass

    def run_dock(self, genes: Optional[List[str]] = None, receptor: Optional[str] = None, receptor_dir: str = "./data/receptors", ligand: Optional[str] = None, ligand_dir: Optional[str] = None, output_dir: str = './docking_results/') -> pl.DataFrame:
        """
        A workflow to perform docking using AutoDock Vina.

        Steps:
            1. Download protein structure from AlphaFold using gene names.
            2. Prepare receptor and ligand files using DockPrep.
            3. Run AutoDock Vina for docking.
            4. Extract and summarize binding affinities.

        Parameters:res
        - receptor: Path to a single receptor PDB file.
        - receptor_dir: Directory containing multiple receptor PDB files.
        - genes: A single gene or a list of gene names to fetch protein structures from AlphaFold.
        - ligand: Path to a single ligand PDB file.
        - ligand_dir: Directory containing multiple ligand PDB files.
        - output_dir: Directory to save docking results.
        """
 
        # Workflow 1: Docking - Ligand + Gene
        # Step 1: Download protein structure from AlphaFold using gene names
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        if not os.path.exists(receptor_dir):
            os.makedirs(receptor_dir, exist_ok=True)

        alphafold = AlphaFold(pdb_path=receptor_dir)

        if len(genes) == 1:
            receptor = alphafold._from_gene_single(genes[0])
        elif len(genes) > 1:
            alphafold._from_gene_batch(genes)
        else:
            print("No genes provided for AlphaFold retrieval.")
            return None
        
        # Step 2: Prepare ligands and proteins
        dockprep = DockPrep()

        # Prepare single ligand
        prepared_ligand = dockprep.prepare_ligand(
        query=ligand,
        target="./data/prepared_ligands/lig1.pdbqt"
        )

        # Prepare single receptor
        prepared_receptor = dockprep.prepare_receptor(
            query=receptor,
            target_prefix="./data/prepared_receptors/rec1",
            AlphaFold=True,
            box_center=None,
            box_size=None
        )

        # print(prepared_ligand, prepared_receptor)

        # Step 3: Docking
        docker = AutoDock()

        # ligand = prepared_ligand
        # receptor = prepared_receptor

        result = docker.singleLigandSingleReceptor(ligand=prepared_ligand, receptor=prepared_receptor, output_dir=output_dir)

        return result

        

