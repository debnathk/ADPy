import os
import argparse
import polars as pl
from typing import Optional

class DataUtils:
    def __init__(self) -> None:
        """Initialize DataUtils class"""
        pass

def extractBindingAffinity(output_file: str) -> str:
    '''
    Extract binding affinity from output .pdbqt file

    Args:
        output_file: Path to PDBQT file generated after docking analysis
        
    Returns:
        Binding_Affinity: Extracted binding affinity in kcal/mol as string, or None if not found

    Raises:
        FileNotFoundError: If the output file doesn't exist
    '''

    with open(output_file, "r", encoding="utf-8", errors="ignore") as f:
        try:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    try:
                        affinity = line.split()[3]
                        break
                    except (IndexError, ValueError) as e:
                        print(f"Error parsing affinity from line '{line.strip()}': {str(e)}")
                        continue

        except IOError as e:
            print(f"Error reading file {output_file}: {str(e)}")

    return affinity

def trimName(filepath):
    """
    Extract filename without Extension from filepath

    Args:

    Returns:

    Example:
        '/path/to/ligand.pdbqt' -> 'ligand'
    """
    filename = os.path.basename(filepath)
    return filename[:-6]

# def log2csv(log_file):

#     '''
#     Save docking output as csv

#     Args:
#         Input: Log file generated after performing docking analysis
#         Output: CSV file containing ligand, receptor and binding affinity
#     '''

#     # Folder containing your .log files
#     # results_folder = output_dir   # change this to your actual folder

#     # Output CSV file
#     output_csv = "example_docking_results.csv"

#     # output_path = os.path.join(output_dir, output_csv)

#     # Collect results
#     results = []

#     # for filename in os.listdir(results_folder):
#     # if filename.endswith("output.log"):
#     #     filepath = os.path.join(results_folder, filename)
#     #     print(filepath)
#     ligand = receptor = affinity = None

#     with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
#         for line in f:
#             if line.startswith("Receptor:"):
#                 receptor = os.path.basename(line.split(":")[1].strip())
#                 # print(receptor)
#             elif line.startswith("Ligand:"):
#                 ligand = os.path.basename(line.split(":")[1].strip())
#                 # print(ligand)
#             elif line.strip().startswith("1"):
#                 parts = line.split()
#                 if len(parts) >= 2:
#                     affinity = parts[1]
#                     break  # stop after first mode (top affinity)
#                 # print(affinity)

#         # if ligand and receptor and affinity:
#         # results.append([ligand, receptor, affinity])

#         df = pd.DataFrame({'Ligand': ligand, 'Receptor': receptor, 'Binding_Affinity(kcal/mol)': affinity}, index=[0])

#     # # Write to CSV
#     # with open(output_path, "w", newline="") as csvfile:
#     #     writer = csv.writer(csvfile)
#     #     writer.writerow(["Ligand", "Receptor", "Affinity (kcal/mol)"])
#     #     writer.writerows(results)

#     # print(f"Results saved to {output_path}")
#     return df

def getAlphaFold(gene: Optional[str] = None, uniprot: Optional[str] = None):
    if gene:
        pass
    if uniprot:
        pass

if __name__ == "__main__":
    pass