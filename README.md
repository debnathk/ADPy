# ADPy: A High-Throughput Tool for Molecular Docking using AutoDock Vina

## Clone Repo

```
git clone https://github.com/debnathk/ADPy.git
cd ADPy
```

## Conda Installation

```
conda env create -f environment.yml
conda activate adpy
```

## Run Docking

### Example

#### Type I: One Ligand - One Receptor Docking

```
python main.py \
    --ligand ./example/prepared_ligands/example_prepared_ligand.pdbqt \
    --receptor ./example/prepared_receptors/example_prepared_receptor.pdbqt \
    --output-dir ./outputs
```

### Ligands and Proteins from User

**Note: Both the ligand and receptor files are needed to be prepared and must be in .pdbqt format.**

#### Type I: One Ligand - One Receptor Docking

```
python main.py --ligand path/to/ligand \
		--receptor path/to/receptor \
		--output-dir path/to/output
```
