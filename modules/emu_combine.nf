/*
========================================================================================
    MODULE : EMU_COMBINE — Tables multi-échantillons par base de données
========================================================================================
    "emu combine-outputs" fusionne les tables de tous les échantillons d'une
    même base en une seule table (une colonne par échantillon), pour chaque
    rang demandé (--emu_combine_ranks), en abondances et en comptages.

    Sorties :
    - {base}/                                   → tables par échantillon (lu par le rapport)
    - results/emu/{base}/emu-combined-{rang}.tsv        → abondances relatives
    - results/emu/{base}/emu-combined-{rang}-counts.tsv → comptages estimés
========================================================================================
*/

process EMU_COMBINE {

    tag "${db_name}"
    label 'process_mono'

    // LEÇON : "pattern" limite la publication aux tables combinées ; les
    // tables par échantillon sont déjà publiées par EMU_ABUNDANCE.
    publishDir "${params.resultsdir}/emu/${db_name}", mode: 'copy', pattern: "emu-combined-*.tsv"

    input:
    tuple val(db_name), path(tables, stageAs: 'tables/*')
    val ranks

    output:
    path "${db_name}",              emit: db_dir
    path "emu-combined-*.tsv",      emit: combined

    script:
    """
    set -euo pipefail

    mkdir "${db_name}"
    cp -L tables/*_rel-abundance.tsv "${db_name}/"

    for rank in ${ranks}; do
        emu combine-outputs "${db_name}" "\$rank"
        emu combine-outputs "${db_name}" "\$rank" --counts
    done

    # Les tables combinées sortent du dossier : il ne garde que les tables
    # par échantillon, lues par le rapport.
    mv "${db_name}"/emu-combined-*.tsv .
    """

    stub:
    """
    mkdir "${db_name}"
    cp -L tables/*_rel-abundance.tsv "${db_name}/"
    for rank in ${ranks}; do
        touch "emu-combined-\$rank.tsv" "emu-combined-\$rank-counts.tsv"
    done
    """
}
