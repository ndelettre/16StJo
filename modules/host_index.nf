/*
========================================================================================
    MODULE : HOST_INDEX — Indexation du génome de l'hôte (minimap2)
========================================================================================
    Construit l'index minimap2 (.mmi) du génome de l'hôte UNE seule fois par
    run, au début (en parallèle de Chopper). Chaque échantillon est ensuite
    aligné sur cet index (HOST_FILTER) sans avoir à réindexer le génome.

    L'index (~7 Go pour GRCh38) reste dans le dossier de travail de Nextflow
    pendant le run, puis est supprimé par HOST_INDEX_CLEANUP dès que tous les
    échantillons sont alignés : rien n'est conservé pour le run suivant.

    Ce module n'est pas utilisé si --host_reference est déjà un index .mmi.

    Sorties :
    - host.mmi → index minimap2 (preset map-ont)
========================================================================================
*/

process HOST_INDEX {

    label 'process_host_index'

    input:
    path host_fasta     // FASTA (.fa/.fasta/.fna, éventuellement .gz)

    output:
    path "host.mmi", emit: index

    script:
    """
    minimap2 -x map-ont -t ${task.cpus} -d host.mmi "${host_fasta}"
    """
}
