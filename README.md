# ADPy: A Python Package for High-Throughput Molecular Docking using AutoDock Vina

## Clone Repo

```
git clone https://github.com/debnathk/ADPy.git
cd ADPy
```

## Setup Virtual Environment (Conda)

```
conda env create -f environment.yml
conda activate adpy
```

## Run Docking

### Example

#### Type I: One Ligand - One Receptor Docking

```
python main.py \
    --ligand ./examples/prepared_ligands/example_prepared_ligand1.pdbqt \
    --receptor ./examples/prepared_receptors/example_prepared_receptor.pdbqt \
    --output-dir ./outputs
```

#### Type II: Multi Ligand - One Receptor Docking

```
python main.py \
    --ligand-dir ./examples/prepared_ligands/ \
    --receptor ./examples/prepared_receptors/example_prepared_receptor.pdbqt \
    --output-dir ./outputs
```

### Docking for User Data

> Please pass appropriate arguments according to your file destinations

**Note: Both the ligand and receptor files are needed to be prepared and must be in .pdbqt format.**
