/*
========================================================================================
    MODULE : CHOPPER — Filtrage qualité et longueur des reads
========================================================================================
    chopper (successeur de NanoFilt) conserve les reads :
      - de qualité moyenne ≥ --min_quality
      - de longueur comprise entre --min_length et --max_length
    Pour le 16S complet (~1 500 pb), la fenêtre de longueur élimine les
    chimères, les fragments et les dimères d'amorces.

    Entrée : un ou plusieurs FASTQ par échantillon (gz ou non) → concaténés
    Sorties :
    - {sample_id}_filtered.fastq.gz → reads filtrés
    - n_reads                       → nombre de reads conservés (pour exclure
                                      les échantillons vides, ex : témoins)
    - {sample_id}_chopper.log       → bilan du filtrage
========================================================================================
*/

process CHOPPER {

    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.resultsdir}/chopper", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    // LEÇON : env() récupère la valeur d'une variable bash définie dans le
    // script. Elle sert dans main.nf à écarter les échantillons trop pauvres.
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), env(N_READS), emit: reads
    path "${sample_id}_chopper.log",                                          emit: log

    script:
    // LEÇON : "zcat -f" lit indifféremment les fichiers compressés ou non,
    // et concatène tous les fichiers d'un dossier barcodeXX en un seul flux.
    """
    set -euo pipefail

    zcat -f ${fastq} \\
        | chopper \\
            --quality ${params.min_quality} \\
            --minlength ${params.min_length} \\
            --maxlength ${params.max_length} \\
            --threads ${task.cpus} \\
            ${params.chopper_args} \\
            2> chopper_stderr.txt \\
        | gzip > "${sample_id}_filtered.fastq.gz"

    N_RAW=\$(zcat -f ${fastq} | awk 'END { print int(NR / 4) }')
    N_READS=\$(zcat "${sample_id}_filtered.fastq.gz" | awk 'END { print int(NR / 4) }')

    {
        echo "Échantillon      : ${sample_id}"
        echo "Fichiers         : ${fastq}"
        echo "Reads bruts      : \$N_RAW"
        echo "Reads conservés  : \$N_READS"
        echo "Qualité min.     : ${params.min_quality}"
        echo "Longueur         : ${params.min_length} - ${params.max_length} pb"
        echo "--- chopper ---"
        cat chopper_stderr.txt
    } > "${sample_id}_chopper.log"

    if [ "\$N_READS" -eq 0 ]; then
        echo "ATTENTION : aucun read n'a survécu au filtrage pour ${sample_id}" >&2
    fi
    """

    stub:
    """
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > "${sample_id}_filtered.fastq.gz"
    N_READS=${sample_id.contains('empty') ? 0 : 1000}
    touch "${sample_id}_chopper.log"
    """
}
