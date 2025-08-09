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

## Usage

### Run Docking

#### Type I: One Ligand - One Receptor Docking

```
python main.py --ligand path/to/ligand --receptor path/to/receptor --output-dir path/to/output
```
