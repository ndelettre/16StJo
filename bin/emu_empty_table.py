#!/usr/bin/env python3
"""
emu_empty_table.py
------------------
Écrit une table d'abondance au format EMU (rel-abundance.tsv, avec la colonne
"estimated counts") pour un échantillon sans aucun read classé : tous ses
reads sont comptés comme "unmapped".

Utilisé quand EMU s'arrête avec "0 reads mapped" (témoin vide, barcode sans
reads...), pour que l'échantillon reste présent dans les tables et le rapport.

Usage :
    emu_empty_table.py --taxonomy <db>/taxonomy.tsv --unmapped 42 \
        --output sample_rel-abundance.tsv
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxonomy", required=True, help="taxonomy.tsv de la base EMU")
    parser.add_argument("--unmapped", required=True, type=int, help="nombre de reads non classés")
    parser.add_argument("--output", required=True, help="fichier TSV de sortie")
    args = parser.parse_args()

    # Mêmes colonnes qu'EMU : tax_id, abundance, rangs de la base, estimated counts
    with open(args.taxonomy, encoding="utf-8") as handle:
        ranks = handle.readline().rstrip("\n").split("\t")[1:]

    header = ["tax_id", "abundance"] + ranks + ["estimated counts"]
    empty_ranks = [""] * len(ranks)
    rows = [
        ["unmapped", "0.0"] + empty_ranks + [f"{float(args.unmapped)}"],
        ["mapped_filtered", "0.0"] + empty_ranks + ["0.0"],
        ["mapped_unclassified", "0.0"] + empty_ranks + ["0.0"],
    ]

    with open(args.output, "w", encoding="utf-8") as out:
        for row in [header] + rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
