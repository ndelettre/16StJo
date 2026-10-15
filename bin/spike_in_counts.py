#!/usr/bin/env python3
"""
spike_in_counts.py
------------------
Compte, pour chaque échantillon, les reads attribués par EMU aux souches du
contrôle interne ZymoBIOMICS (Spike-in Control I ou II).

Pour chaque échantillon et chaque souche attendue :
  - estimated_reads : nombre de reads estimé par EMU
  - pct_classified  : % des reads classés (bactéries identifiées)
  - pct_total       : % de tous les reads analysés par EMU (classés + non classés)

Une ligne "<genre> (autres espèces)" est ajoutée quand des reads sont
attribués à une autre espèce d'un genre du contrôle (ex. Allobacillus salarius) :
signe qu'une partie du contrôle a pu être mal attribuée.

Usage :
    spike_in_counts.py --control "Spike-in I" \
        --species "Imtechella halotolerans,Allobacillus halotolerans" \
        --output spike_in_counts.tsv  sample1_rel-abundance.tsv ...
"""

import argparse
import os

import pandas as pd

NOT_CLASSIFIED = {"unmapped", "mapped_filtered", "mapped_unclassified"}


def sample_name(path):
    return os.path.basename(path).replace("_rel-abundance.tsv", "")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--control", required=True, help="nom du contrôle interne")
    parser.add_argument("--species", required=True, help="espèces attendues, séparées par des virgules")
    parser.add_argument("--output", required=True, help="fichier TSV de sortie")
    parser.add_argument("tables", nargs="+", help="tables EMU {sample}_rel-abundance.tsv")
    args = parser.parse_args()

    expected = [s.strip() for s in args.species.split(",") if s.strip()]
    expected_genera = {s.split()[0] for s in expected}

    rows = []
    for path in sorted(args.tables, key=sample_name):
        sample = sample_name(path)
        df = pd.read_csv(path, sep="\t", dtype={"tax_id": str})
        counts = pd.to_numeric(df["estimated counts"], errors="coerce").fillna(0)
        classified = ~df["tax_id"].isin(NOT_CLASSIFIED)
        total_classified = counts[classified].sum()
        total = counts.sum()

        def pct(value, denominator):
            return round(100 * value / denominator, 2) if denominator > 0 else 0.0

        for species in expected:
            reads = counts[classified & (df["species"] == species)].sum()
            rows.append([sample, args.control, species, "oui", round(reads, 1),
                         pct(reads, total_classified), pct(reads, total)])

        for genus in sorted(expected_genera):
            others = classified & (df["genus"] == genus) & ~df["species"].isin(expected)
            reads = counts[others].sum()
            if reads > 0:
                rows.append([sample, args.control, f"{genus} (autres espèces)", "non",
                             round(reads, 1), pct(reads, total_classified), pct(reads, total)])

    out = pd.DataFrame(rows, columns=["sample", "control", "species", "expected",
                                      "estimated_reads", "pct_classified", "pct_total"])
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
