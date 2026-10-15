/*
========================================================================================
    MODULE : HOST_INDEX_CLEANUP — Suppression de l'index du génome de l'hôte
========================================================================================
    Supprime l'index .mmi construit par HOST_INDEX (~7 Go pour GRCh38) dès que
    TOUS les échantillons ont été alignés sur l'hôte, sans attendre la fin
    d'EMU. N'est jamais appelé sur un .mmi fourni par l'utilisateur.

    Conséquence pour "-resume" : l'index étant supprimé, il sera reconstruit
    si le run est relancé, et les alignements sur l'hôte seront refaits.
========================================================================================
*/

process HOST_INDEX_CLEANUP {

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : un process "exec:" exécute du code Groovy directement dans
    // Nextflow (pas de conteneur, pas de script shell). L'index est reçu en
    // "val" et non en "path" : on obtient ainsi son chemin réel, dans le
    // dossier de travail de HOST_INDEX, au lieu d'une copie locale.
    //
    // LEÇON : le second input (`all_done`) ne sert qu'à attendre : il reçoit la
    // liste de tous les résultats de HOST_FILTER (via .collect()), qui n'est
    // émise qu'une fois le dernier échantillon aligné.
    // ─────────────────────────────────────────────────────────────────────────
    input:
    val index
    val all_done

    exec:
    def index_file = index.toRealPath()
    def size_gb = String.format('%.1f', index_file.size() / 1e9)
    index_file.delete()
    log.info "Index du génome hôte supprimé (${size_gb} Go) : ${index_file}"
}
