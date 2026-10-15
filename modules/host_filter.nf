/*
========================================================================================
    MODULE : HOST_FILTER — Comptage et retrait des reads de l'hôte
========================================================================================
    Pour un échantillon : compte les reads identifiés comme "hôte" par HOST_MAP,
    puis écrit un FASTQ sans ces reads (c'est lui qui part dans EMU).

    Sorties :
    - {sample_id}_nohost.fastq.gz → reads sans l'hôte (entrée d'EMU)
    - {sample_id}_host.tsv        → sample | reads_filtered | host_reads | host_pct | reads_nohost
========================================================================================
*/

process HOST_FILTER {

    tag "${sample_id}"
    label 'process_mono'

    publishDir "${params.resultsdir}/host_removed", mode: 'copy', pattern: "*_nohost.fastq.gz"

    input:
    tuple val(sample_id), path(fastq)
    path host_reads     // host_reads.tsv produit par HOST_MAP (commun à tous les échantillons)

    output:
    tuple val(sample_id), path("${sample_id}_nohost.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_host.tsv"),        emit: stats

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : la liste des reads hôte est chargée dans un tableau awk (h[])
    // dans le bloc BEGIN avec getline. On n'utilise pas l'astuce classique
    // "NR == FNR" car elle se trompe quand le premier fichier est vide
    // (échantillon sans aucun read hôte).
    // ─────────────────────────────────────────────────────────────────────────
    script:
    """
    set -euo pipefail

    awk -F'\\t' -v s="${sample_id}" '\$1 == s { print \$2 }' "${host_reads}" > ids.txt

    gzip -cd "${fastq}" \\
    | awk -v ids=ids.txt -v sample="${sample_id}" -v stats="${sample_id}_host.tsv" '
        BEGIN { while ((getline id < ids) > 0) h[id] = 1 }
        NR % 4 == 1 { total++; keep = !(substr(\$1, 2) in h); if (!keep) host++ }
        keep { print }
        END {
            pct = (total > 0) ? 100 * host / total : 0
            print "sample\\treads_filtered\\thost_reads\\thost_pct\\treads_nohost" > stats
            printf "%s\\t%d\\t%d\\t%.2f\\t%d\\n", sample, total, host, pct, total - host > stats
        }
      ' \\
    | gzip > "${sample_id}_nohost.fastq.gz"

    cat "${sample_id}_host.tsv" >&2
    """
}
