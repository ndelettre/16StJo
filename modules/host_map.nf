/*
========================================================================================
    MODULE : HOST_MAP — Alignement de tous les reads contre le génome de l'hôte
========================================================================================
    minimap2 aligne les reads filtrés de TOUS les échantillons contre le génome
    de l'hôte (GRCh38 par défaut) en un seul passage.

    Pourquoi un seul passage pour tous les échantillons ?
    Charger ou indexer le génome humain prend plusieurs minutes et ~10-15 Go
    de RAM. Le faire une fois par échantillon (jusqu'à 94) serait très lent.
    Ici l'index est construit (ou chargé, si un .mmi est fourni) UNE seule
    fois, en mémoire : aucun fichier d'index n'est écrit sur le disque.

    Pour retrouver l'échantillon d'origine de chaque read, son nom est préfixé
    par le numéro de l'échantillon ("3__<read_id>") avant l'alignement.

    Un read est compté comme "hôte" si son alignement couvre au moins 50 %
    de sa longueur : on évite ainsi de compter un read bactérien qui n'aurait
    qu'un court alignement fortuit sur le génome humain.

    Sorties :
    - host_reads.tsv → une ligne par read hôte : <sample_id> <TAB> <read_id>
========================================================================================
*/

process HOST_MAP {

    label 'process_host'

    input:
    path reads          // tous les {sample_id}_filtered.fastq.gz
    path host_reference // FASTA (.fa/.fasta/.fna, éventuellement .gz) ou index .mmi

    output:
    path "host_reads.tsv", emit: hits

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : minimap2 lit les reads sur son entrée standard grâce à "-".
    // La boucle "for" écrit les reads de chaque échantillon les uns après les
    // autres dans ce flux, après avoir renommé chaque read avec awk
    // (NR % 4 == 1 → 1re ligne de chaque bloc de 4 lignes FASTQ = l'en-tête).
    //
    // Colonnes PAF utilisées : $1 = nom du read, $2 = longueur du read,
    // $3/$4 = début/fin de la partie alignée du read.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    """
    set -euo pipefail

    i=0
    : > samples_index.tsv
    for fq in *_filtered.fastq.gz; do
        sample_id="\${fq%_filtered.fastq.gz}"
        printf '%s\\t%s\\n' "\$i" "\$sample_id" >> samples_index.tsv
        gzip -cd "\$fq" | awk -v prefix="\$i" 'NR % 4 == 1 { \$0 = "@" prefix "__" substr(\$1, 2) } { print }'
        i=\$((i + 1))
    done \\
    | minimap2 -x map-ont --secondary=no -t ${task.cpus} "${host_reference}" - \\
    > host.paf

    # Reads hôte : alignement couvrant ≥ 50 % du read → "<sample_id>\\t<read_id>"
    awk -F'\\t' '(\$4 - \$3) >= 0.5 * \$2 { print \$1 }' host.paf \\
    | sort -u \\
    | awk -F'\\t' -v index_file=samples_index.tsv '
        BEGIN { while ((getline line < index_file) > 0) { split(line, f, "\\t"); name[f[1]] = f[2] } }
        { pos = index(\$0, "__"); print name[substr(\$0, 1, pos - 1)] "\\t" substr(\$0, pos + 2) }
      ' > host_reads.tsv

    rm -f host.paf
    echo "Reads hôte détectés (tous échantillons) : \$(wc -l < host_reads.tsv)" >&2
    """
}
