"""
Prepare the receptor structure for docking.

Downloads PDB 1IEP (Abl kinase domain in complex with imatinib), strips
waters/heteroatoms and the co-crystallized ligand, and converts the
resulting protein-only structure to PDBQT (adding polar hydrogens and
Gasteiger charges) using OpenBabel.

Also records the docking search box, centered on the original imatinib
binding pocket (taken from the 1IEP ligand coordinates), so dock.py can
reuse it.
"""

import json
import os
import subprocess

import requests

PDB_ID = "1IEP"
STRUCT_DIR = "structures"
RAW_PDB = os.path.join(STRUCT_DIR, f"{PDB_ID}.pdb")
PROTEIN_PDB = os.path.join(STRUCT_DIR, "receptor_protein.pdb")
RECEPTOR_PDBQT = os.path.join(STRUCT_DIR, "receptor.pdbqt")
BOX_PATH = os.path.join(STRUCT_DIR, "box.json")

# Search box centered on the imatinib binding site in 1IEP (from published
# docking tutorials using this exact benchmark structure), sized generously
# to cover the ATP-binding pocket.
DOCKING_BOX = {
    "center_x": 15.190,
    "center_y": 53.903,
    "center_z": 16.917,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,
}


def download_pdb():
    os.makedirs(STRUCT_DIR, exist_ok=True)
    if os.path.exists(RAW_PDB):
        print(f"  {RAW_PDB} already downloaded")
        return
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    with open(RAW_PDB, "w") as f:
        f.write(resp.text)
    print(f"  downloaded {PDB_ID} -> {RAW_PDB}")


def strip_to_protein_only():
    """Keep only ATOM records (protein backbone/sidechains), drop waters,
    ions, and the co-crystallized ligand (HETATM records)."""
    with open(RAW_PDB) as f:
        lines = f.readlines()

    kept = [l for l in lines if l.startswith("ATOM") or l.startswith("TER") or l.startswith("END")]

    with open(PROTEIN_PDB, "w") as f:
        f.writelines(kept)
    print(f"  stripped to protein-only atoms -> {PROTEIN_PDB} ({len(kept)} lines)")


def convert_to_pdbqt():
    """Use OpenBabel to add polar hydrogens, assign Gasteiger charges, and
    write a receptor PDBQT suitable for AutoDock Vina."""
    cmd = [
        "obabel",
        PROTEIN_PDB,
        "-O", RECEPTOR_PDBQT,
        "-xr",       # rigid receptor
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"obabel failed:\n{result.stderr}")
    print(f"  converted to PDBQT -> {RECEPTOR_PDBQT}")


def write_box():
    with open(BOX_PATH, "w") as f:
        json.dump(DOCKING_BOX, f, indent=2)
    print(f"  wrote docking search box -> {BOX_PATH}")


def main():
    print(f"Preparing receptor: PDB {PDB_ID} (Abl kinase domain / imatinib complex)")
    download_pdb()
    strip_to_protein_only()
    convert_to_pdbqt()
    write_box()
    print("\nReceptor ready.")


if __name__ == "__main__":
    main()
