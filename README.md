# Molecular Docking Pipeline for Drug Discovery

A small end-to-end cheminformatics pipeline: pull candidate small-molecule
structures, clean/featurize them with RDKit, dock them against a target
protein with AutoDock Vina, and rank candidates by predicted binding
affinity.

This mirrors an early-stage virtual screening step in a real drug discovery
workflow — the same task that sits behind terms like "molecular modeling,"
"docking," and "cheminformatics" in computational drug discovery job
postings.

## Target system

- **Receptor:** Abl tyrosine kinase, PDB [1IEP](https://www.rcsb.org/structure/1IEP)
  (kinase domain bound to imatinib/Gleevec) — a well-characterized,
  commonly-used docking benchmark target.
- **Ligand library:** a curated set of small molecules (including imatinib
  itself as a positive control) pulled live from PubChem by name/CID, plus
  a fallback offline CSV so the pipeline runs without network access.

## Pipeline

```
fetch_ligands.py     Pull candidate compounds from PubChem -> data/ligands.csv (SMILES)
prepare_ligands.py    RDKit: parse SMILES, embed 3D, minimize, compute Lipinski
                       descriptors, filter, write PDBQT per ligand (via Meeko)
prepare_receptor.py   Download 1IEP, strip waters/heteroatoms, protonate,
                       convert to PDBQT (via OpenBabel), define search box
dock.py                Run AutoDock Vina for every prepared ligand against the
                       receptor, collect binding-affinity scores
rank_results.py        Merge docking scores + RDKit descriptors, output a
                       ranked leaderboard (results/ranked_candidates.csv)
```

## Setup

```bash
conda env create -f environment.yml
conda activate docking
```

## Run

```bash
python fetch_ligands.py
python prepare_receptor.py
python prepare_ligands.py
python dock.py
python rank_results.py
```

Output: `results/ranked_candidates.csv` — candidate compounds ranked by
predicted binding affinity (kcal/mol, more negative = stronger predicted
binding), alongside Lipinski Rule-of-Five descriptors (MW, LogP, HBD, HBA,
TPSA) for quick druglikeness triage.

## Why these choices

- **RDKit** for all cheminformatics (SMILES parsing, 3D embedding, Lipinski
  descriptors) — the standard open-source cheminformatics toolkit.
- **AutoDock Vina** for docking — fast, widely used, and has first-class
  Python bindings (`vina` package), so the whole pipeline stays in Python.
- **Meeko** for ligand PDBQT preparation (the modern AutoDockTools
  replacement) and **OpenBabel** for receptor protonation/conversion.
- **1IEP** was chosen as the target because it's a validated, literature-
  documented docking benchmark (Abl kinase + imatinib), which lets the
  pipeline's output be sanity-checked against a known active compound.

## Status

Built as a portfolio project to demonstrate computational drug-discovery
skills (molecular modeling, docking, cheminformatics, chemical data
pipelines) for internship applications in modeling & informatics /
computational chemistry.
