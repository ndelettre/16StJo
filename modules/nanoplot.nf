/*
========================================================================================
    MODULE : NANOPLOT — Contrôle qualité des reads
========================================================================================
    NanoPlot produit un rapport HTML interactif (distribution des longueurs,
    des qualités, longueur vs qualité...) et un fichier de statistiques.
    Il est lancé deux fois (alias dans main.nf) :
      - stage "raw"      → reads bruts, avant chopper
      - stage "filtered" → reads conservés après chopper

    Sorties :
    - {sample_id}_{stage}/                  → rapport NanoPlot complet
    - {sample_id}_{stage}_NanoStats.txt     → statistiques (lues par le rapport)
========================================================================================
*/

process NANOPLOT {

    tag "${sample_id} (${stage})"
    label 'process_low'

    publishDir "${params.resultsdir}/nanoplot", mode: 'copy', pattern: "${sample_id}_${stage}"

    input:
    tuple val(sample_id), val(stage), path(fastq)

    output:
    path "${sample_id}_${stage}",                  emit: report
    path "${sample_id}_${stage}_NanoStats.txt",    emit: stats

    script:
    // LEÇON : --no_static désactive les images PNG (qui nécessitent un
    // navigateur dans le conteneur) ; le rapport HTML reste interactif.
    // --tsv_stats écrit NanoStats.txt en tableau "métrique<TAB>valeur".
    // MPLCONFIGDIR : le conteneur tourne avec l'UID de l'utilisateur
    // (docker -u), qui n'a pas de HOME inscriptible.
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export HOME=\$PWD

    NanoPlot \\
        --fastq ${fastq} \\
        --threads ${task.cpus} \\
        --outdir "${sample_id}_${stage}" \\
        --title "${sample_id} (${stage})" \\
        --no_static \\
        --tsv_stats

    cp "${sample_id}_${stage}/NanoStats.txt" "${sample_id}_${stage}_NanoStats.txt"
    """

    stub:
    """
    mkdir -p "${sample_id}_${stage}"
    printf 'Metrics\\tdataset\\nnumber_of_reads\\t1000.0\\nnumber_of_bases\\t1500000.0\\nmedian_read_length\\t1500.0\\nmean_read_length\\t1490.0\\nmedian_qual\\t14.0\\nmean_qual\\t13.5\\nn50\\t1510.0\\n' \\
        > "${sample_id}_${stage}/NanoStats.txt"
    cp "${sample_id}_${stage}/NanoStats.txt" "${sample_id}_${stage}_NanoStats.txt"
    """
}
