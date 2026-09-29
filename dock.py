"""
Dock every prepared ligand against the prepared receptor using AutoDock
Vina's Python bindings, within the search box defined in
structures/box.json.

Output: data/docking_scores.csv (name, best predicted binding affinity in
kcal/mol; more negative = stronger predicted binding).
"""

import csv
import glob
import json
import os

from vina import Vina

RECEPTOR_PDBQT = os.path.join("structures", "receptor.pdbqt")
BOX_PATH = os.path.join("structures", "box.json")
LIGAND_DIR = os.path.join("structures", "ligands")
SCORES_CSV = os.path.join("data", "docking_scores.csv")


def dock_ligand(v: Vina, ligand_path: str) -> float:
    v.set_ligand_from_file(ligand_path)
    v.dock(exhaustiveness=8, n_poses=9)
    energies = v.energies(n_poses=1)
    return float(energies[0][0])  # best pose binding affinity, kcal/mol


def main():
    with open(BOX_PATH) as f:
        box = json.load(f)

    v = Vina(sf_name="vina")
    v.set_receptor(RECEPTOR_PDBQT)
    v.compute_vina_maps(
        center=[box["center_x"], box["center_y"], box["center_z"]],
        box_size=[box["size_x"], box["size_y"], box["size_z"]],
    )

    ligand_paths = sorted(glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt")))
    if not ligand_paths:
        raise SystemExit("No prepared ligands found — run prepare_ligands.py first.")

    print(f"Docking {len(ligand_paths)} ligands against {RECEPTOR_PDBQT}...")
    rows = []
    for path in ligand_paths:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            affinity = dock_ligand(v, path)
            rows.append({"name": name, "binding_affinity_kcal_mol": affinity})
            print(f"  {name}: {affinity:.2f} kcal/mol")
        except Exception as exc:
            print(f"  WARNING: docking failed for {name}: {exc}")

    with open(SCORES_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "binding_affinity_kcal_mol"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote docking scores for {len(rows)} ligands -> {SCORES_CSV}")


if __name__ == "__main__":
    main()
