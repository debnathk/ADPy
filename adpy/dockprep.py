import os
import csv
from typing import List, Tuple, Optional, Union
import subprocess
import requests
import sys

from vina import Vina
from .utils import trimName, extractBindingAffinity
import polars as pl

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