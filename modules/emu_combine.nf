/*
========================================================================================
    MODULE : EMU_COMBINE — Tables multi-échantillons par rang taxonomique
========================================================================================
    `emu combine-outputs` fusionne les tables de tous les échantillons en une
    table par rang taxonomique (une colonne par échantillon), en abondance
    relative et en nombre de reads estimé.

    Sorties (pour rang = species, genus, family) :
    - emu-combined-{rang}.tsv         → abondances relatives
    - emu-combined-{rang}-counts.tsv  → nombres de reads estimés

    Note : dans ces tables, la ligne sans nom de taxon regroupe les reads non
    classés (unmapped + mapped_filtered + mapped_unclassified).
========================================================================================
*/

process EMU_COMBINE {

    label 'process_low'

    publishDir "${params.resultsdir}/emu", mode: 'copy'

    input:
    path abundance_tables   // tous les {sample_id}_rel-abundance.tsv

    output:
    path "emu-combined-*.tsv", emit: tables

    script:
    """
    set -euo pipefail

    for rank in species genus family; do
        emu combine-outputs . "\$rank"
        emu combine-outputs . "\$rank" --counts
    done
    """
}
