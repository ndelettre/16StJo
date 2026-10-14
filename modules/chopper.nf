/*
========================================================================================
    MODULE : CHOPPER — Filtrage des reads par longueur et qualité
========================================================================================
    Chopper (successeur de NanoFilt) garde les reads dont :
      - la longueur est comprise entre min_length et max_length
      - le Q-score moyen est ≥ min_quality

    Pour du 16S full-length (~1 500 bp), la fenêtre 1 200 - 1 800 bp écarte
    les fragments tronqués et les chimères/concatémères trop longs.

    Sorties :
    - {sample_id}_filtered.fastq.gz → reads conservés
    - {sample_id}_chopper.log       → "Kept X reads out of Y reads"
========================================================================================
*/

process CHOPPER {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.resultsdir}/chopper", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_chopper.log"),       emit: log

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : Chopper lit directement les .fastq et .fastq.gz (option -i)
    // et écrit les reads conservés sur stdout ; on les compresse avec gzip.
    // Son bilan ("Kept X reads out of Y reads") part sur stderr : on le
    // redirige dans le fichier log.
    //
    // Un FASTQ d'entrée vide (barcode sans reads) n'est pas une erreur : on
    // produit un FASTQ vide pour que l'échantillon reste visible dans le
    // rapport avec 0 read, au lieu d'arrêter tout le run.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    """
    set -euo pipefail

    if [[ "${fastq}" == *.gz ]]; then
        n_bytes=\$( (gzip -cd "${fastq}" 2>/dev/null || true) | head -c 1 | wc -c)
    else
        n_bytes=\$(head -c 1 "${fastq}" | wc -c)
    fi

    if [ "\$n_bytes" -eq 0 ]; then
        echo "Kept 0 reads out of 0 reads" > "${sample_id}_chopper.log"
        printf '' | gzip > "${sample_id}_filtered.fastq.gz"
    else
        chopper \\
            --quality ${params.min_quality} \\
            --minlength ${params.min_length} \\
            --maxlength ${params.max_length} \\
            --threads ${task.cpus} \\
            --input "${fastq}" \\
            2> "${sample_id}_chopper.log" \\
            | gzip > "${sample_id}_filtered.fastq.gz"
    fi

    cat "${sample_id}_chopper.log" >&2
    """
}
