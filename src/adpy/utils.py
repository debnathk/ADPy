import os
import csv
import argparse
import pandas as pd

def log2csv(log_file):

    # Folder containing your .log files
    # results_folder = output_dir   # change this to your actual folder

    # Output CSV file
    output_csv = "example_docking_results.csv"

    # output_path = os.path.join(output_dir, output_csv)

    # Collect results
    results = []

    # for filename in os.listdir(results_folder):
    # if filename.endswith("output.log"):
    #     filepath = os.path.join(results_folder, filename)
    #     print(filepath)
    ligand = receptor = affinity = None

    with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("Receptor:"):
                receptor = os.path.basename(line.split(":")[1].strip())
                # print(receptor)
            elif line.startswith("Ligand:"):
                ligand = os.path.basename(line.split(":")[1].strip())
                # print(ligand)
            elif line.strip().startswith("1"):
                parts = line.split()
                if len(parts) >= 2:
                    affinity = parts[1]
                    break  # stop after first mode (top affinity)
                # print(affinity)

        # if ligand and receptor and affinity:
        # results.append([ligand, receptor, affinity])

        df = pd.DataFrame({'Ligand': ligand, 'Receptor': receptor, 'Binding_Affinity(kcal/mol)': affinity}, index=[0])

    # # Write to CSV
    # with open(output_path, "w", newline="") as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(["Ligand", "Receptor", "Affinity (kcal/mol)"])
    #     writer.writerows(results)

    # print(f"Results saved to {output_path}")
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--log_file', required=True, help='Set log file path')
    args = parser.parse_args()

    log2csv(log_file=args.log_file)