"""
Prepare candidate ligands for docking.

For each candidate in data/ligands.csv:
  1. Parse the SMILES with RDKit and add explicit hydrogens.
  2. Embed a 3D conformer and run a quick MMFF minimization.
  3. Compute Lipinski Rule-of-Five descriptors (MW, LogP, HBD, HBA, TPSA)
     for a fast druglikeness triage.
  4. Write a PDBQT file (via Meeko) for AutoDock Vina.

Outputs:
  data/ligand_descriptors.csv   RDKit descriptors per candidate
  structures/ligands/<name>.pdbqt
"""

import csv
import os

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from meeko import MoleculePreparation, PDBQTWriterLegacy

LIGANDS_CSV = os.path.join("data", "ligands.csv")
DESCRIPTORS_CSV = os.path.join("data", "ligand_descriptors.csv")
LIGAND_PDBQT_DIR = os.path.join("structures", "ligands")


def compute_descriptors(mol) -> dict:
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
    }


def lipinski_violations(desc: dict) -> int:
    violations = 0
    if desc["mw"] > 500:
        violations += 1
    if desc["logp"] > 5:
        violations += 1
    if desc["hbd"] > 5:
        violations += 1
    if desc["hba"] > 10:
        violations += 1
    return violations


def prepare_one(name: str, smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  WARNING: RDKit could not parse SMILES for {name}: {smiles}")
        return None

    desc = compute_descriptors(mol)
    desc["name"] = name
    desc["smiles"] = smiles
    desc["lipinski_violations"] = lipinski_violations(desc)

    # 3D embedding + MMFF minimization
    mol_h = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
    if embed_status != 0:
        print(f"  WARNING: 3D embedding failed for {name}")
        return None
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)

    # Meeko: RDKit mol -> docking-ready PDBQT
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol_h)
    pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        print(f"  WARNING: Meeko PDBQT export failed for {name}: {err}")
        return None

    os.makedirs(LIGAND_PDBQT_DIR, exist_ok=True)
    out_path = os.path.join(LIGAND_PDBQT_DIR, f"{name}.pdbqt")
    with open(out_path, "w") as f:
        f.write(pdbqt_string)

    print(f"  prepared {name}: MW={desc['mw']}, LogP={desc['logp']}, "
          f"Lipinski violations={desc['lipinski_violations']}")
    return desc


def main():
    with open(LIGANDS_CSV) as f:
        ligands = list(csv.DictReader(f))

    print(f"Preparing {len(ligands)} ligands...")
    all_desc = []
    for row in ligands:
        desc = prepare_one(row["name"], row["smiles"])
        if desc:
            all_desc.append(desc)

    fieldnames = ["name", "smiles", "mw", "logp", "hbd", "hba", "tpsa",
                  "rotatable_bonds", "lipinski_violations"]
    with open(DESCRIPTORS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_desc)

    print(f"\nWrote descriptors for {len(all_desc)} ligands -> {DESCRIPTORS_CSV}")
    print(f"Wrote PDBQT files -> {LIGAND_PDBQT_DIR}/")


if __name__ == "__main__":
    main()
