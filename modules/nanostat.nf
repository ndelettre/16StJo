/*
========================================================================================
    MODULE : NANOSTAT — Statistiques qualité des reads
========================================================================================
    NanoStat calcule les métriques de qualité d'un FASTQ (nombre de reads,
    longueurs, N50, Q-score moyen/médian...).

    Le module est appelé deux fois dans main.nf (alias NANOSTAT_RAW et
    NANOSTAT_FILTERED) : avant et après Chopper. La valeur "stage"
    ("raw" ou "filtered") distingue les deux passages.

    Sorties :
    - {sample_id}_{stage}_nanostat.tsv → métriques au format TSV (lu par le rapport)
========================================================================================
*/

process NANOSTAT {

    tag "${sample_id} (${stage})"
    label 'process_mono'

    publishDir "${params.resultsdir}/nanostat", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    val stage

    output:
    tuple val(sample_id), path("${sample_id}_${stage}_nanostat.tsv"), emit: stats

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : NanoStat plante ("No reads found in input") sur un FASTQ vide,
    // ce qui arrive pour un barcode sans reads ou si Chopper a tout écarté.
    // Avec 94 échantillons, un seul fichier vide ne doit pas arrêter le run :
    // dans ce cas on écrit nous-mêmes un TSV à 0 read, au même format.
    //
    // LEÇON : `head -c 1` ferme le pipe dès le premier octet lu, ce qui envoie
    // un SIGPIPE à gzip. Le `|| true` évite que `set -o pipefail` transforme
    // ce signal normal en erreur du process.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    def out = "${sample_id}_${stage}_nanostat.tsv"
    """
    set -euo pipefail

    if [[ "${fastq}" == *.gz ]]; then
        n_bytes=\$( (gzip -cd "${fastq}" 2>/dev/null || true) | head -c 1 | wc -c)
    else
        n_bytes=\$(head -c 1 "${fastq}" | wc -c)
    fi

    if [ "\$n_bytes" -eq 0 ]; then
        echo "ATTENTION : aucun read dans ${fastq} (${sample_id}, ${stage})" >&2
        printf 'Metrics\\tdataset\\n' > "${out}"
        for metric in number_of_reads number_of_bases median_read_length \\
                      mean_read_length read_length_stdev n50 mean_qual median_qual; do
            printf '%s\\t0\\n' "\$metric" >> "${out}"
        done
    else
        NanoStat \\
            --fastq "${fastq}" \\
            --tsv \\
            --name "${out}" \\
            --threads ${task.cpus}
    fi
    """
}
