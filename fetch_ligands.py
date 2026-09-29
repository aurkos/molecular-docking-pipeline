"""
Fetch candidate small-molecule structures from PubChem.

Pulls SMILES for a curated candidate list (kinase-inhibitor-like small
molecules, including imatinib as a positive control against the 1IEP
receptor) via the PubChem PUG REST API, and writes them to
data/ligands.csv. Falls back to a bundled offline CSV if PubChem is
unreachable, so the rest of the pipeline never blocks on network access.
"""

import csv
import os
import time

import requests

PUBCHEM_BASE = "https://pubchem.ncbi.gov/rest/pug"

# Curated candidate set: imatinib (known 1IEP binder, positive control)
# plus a handful of other kinase-inhibitor-like small molecules for the
# pipeline to screen and rank.
CANDIDATES = [
    "imatinib",
    "nilotinib",
    "dasatinib",
    "erlotinib",
    "gefitinib",
    "sorafenib",
    "sunitinib",
    "aspirin",       # non-kinase-inhibitor control, expect weak binding
    "caffeine",       # non-kinase-inhibitor control, expect weak binding
]

OUT_PATH = os.path.join("data", "ligands.csv")
FALLBACK_PATH = os.path.join("data", "ligands_offline_fallback.csv")


def fetch_smiles(name: str) -> tuple[str, str] | None:
    """Look up a compound's CID and isomeric SMILES on PubChem by name."""
    url = f"{PUBCHEM_BASE}/compound/name/{name}/property/IsomericSMILES/JSON"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    props = resp.json()["PropertyTable"]["Properties"][0]
    return str(props["CID"]), props["IsomericSMILES"]


def main():
    os.makedirs("data", exist_ok=True)
    rows = []

    for name in CANDIDATES:
        try:
            cid, smiles = fetch_smiles(name)
            rows.append({"name": name, "cid": cid, "smiles": smiles})
            print(f"  fetched {name}: CID {cid}")
            time.sleep(0.2)  # be polite to the public API
        except Exception as exc:
            print(f"  WARNING: could not fetch {name} from PubChem ({exc})")

    if not rows and os.path.exists(FALLBACK_PATH):
        print("PubChem unreachable — using bundled offline fallback data.")
        with open(FALLBACK_PATH) as f:
            rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit("No ligand data available (no network and no fallback file).")

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "cid", "smiles"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} ligands to {OUT_PATH}")


if __name__ == "__main__":
    main()
