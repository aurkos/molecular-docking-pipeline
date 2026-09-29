"""
Merge docking scores with RDKit descriptors into a single ranked
leaderboard of candidates, sorted by predicted binding affinity
(most negative / strongest predicted binder first).
"""

import csv
import os

DESCRIPTORS_CSV = os.path.join("data", "ligand_descriptors.csv")
SCORES_CSV = os.path.join("data", "docking_scores.csv")
OUT_CSV = os.path.join("results", "ranked_candidates.csv")


def main():
    with open(DESCRIPTORS_CSV) as f:
        descriptors = {row["name"]: row for row in csv.DictReader(f)}

    with open(SCORES_CSV) as f:
        scores = {row["name"]: row["binding_affinity_kcal_mol"] for row in csv.DictReader(f)}

    merged = []
    for name, desc in descriptors.items():
        if name not in scores:
            continue
        row = dict(desc)
        row["binding_affinity_kcal_mol"] = float(scores[name])
        merged.append(row)

    merged.sort(key=lambda r: r["binding_affinity_kcal_mol"])

    os.makedirs("results", exist_ok=True)
    fieldnames = ["name", "binding_affinity_kcal_mol", "mw", "logp", "hbd",
                  "hba", "tpsa", "rotatable_bonds", "lipinski_violations", "smiles"]
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged)

    print(f"Ranked {len(merged)} candidates -> {OUT_CSV}\n")
    print(f"{'Rank':<5}{'Name':<15}{'Affinity (kcal/mol)':<22}{'Lipinski violations'}")
    for i, row in enumerate(merged, start=1):
        print(f"{i:<5}{row['name']:<15}{row['binding_affinity_kcal_mol']:<22}{row['lipinski_violations']}")


if __name__ == "__main__":
    main()
