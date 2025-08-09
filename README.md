# ADPy: A High-Throughput Tool for Molecular Docking using AutoDock Vina

# Installation

## Using pip

```
pip install adpy
```

## Using conda

```
conda env create -f environment.yml
conda activate adpy
```

## Usage

### Run Docking

```
from adpy import run_dock

prep_ligand = '../examples/prepared_ligands/example_prepared_ligand.pdbqt'
prep_receptor = '../examples/prepared_proteins/example_prepared_protein.pdbqt'

run_dock(ligand=prep_ligand, receptor=prep_receptor)
```

### Save Output (as .csv)

```
from adpy import pdbqt2csv

pdqt_file = "../examples/results/example_docking.pdbqt"

result = pdbqt2csv(pdbqt_file=pdqt_file)
```
