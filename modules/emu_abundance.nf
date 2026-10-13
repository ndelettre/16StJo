/*
========================================================================================
    MODULE : EMU_ABUNDANCE — Profil taxonomique 16S
========================================================================================
    Emu aligne les reads contre une base 16S (minimap2) puis estime les
    abondances relatives des espèces par un algorithme EM qui tient compte
    du taux d'erreur des reads Nanopore.

    Le process est lancé pour chaque couple (échantillon, base) : les mêmes
    reads sont classés contre chaque base choisie avec --emu_dbs.

    Sorties (results/emu/{base}/) :
    - {sample_id}_rel-abundance.tsv → abondances + comptages estimés par taxon
      + lignes "unmapped", "mapped_filtered", "mapped_unclassified"
    - {sample_id}_rel-abundance-threshold-X.tsv → idem après seuil d'abondance
    - {sample_id}_emu.log           → log Emu (nb reads alignés / non alignés)
========================================================================================
*/

process EMU_ABUNDANCE {

    tag "${sample_id} | ${db_name}"
    label 'process_medium'

    publishDir "${params.resultsdir}/emu/${db_name}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(db_name), path(db)

    output:
    tuple val(db_name), val(sample_id), path("${sample_id}_rel-abundance.tsv"), emit: abundance
    path "${sample_id}_rel-abundance-threshold-*.tsv", optional: true,         emit: threshold
    path "${sample_id}_emu.log",                                               emit: log

    script:
    // LEÇON : --keep-counts ajoute la colonne "estimated counts" (nb de reads
    // estimé par taxon), indispensable pour les tables de comptage et pour
    // calculer le taux de classification dans le rapport.
    """
    set -euo pipefail

    emu abundance \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --min-abundance ${params.emu_min_abundance} \\
        --N ${params.emu_N} \\
        --keep-counts \\
        --output-dir . \\
        --output-basename ${sample_id} \\
        ${params.emu_args} \\
        ${reads} \\
        > "${sample_id}_emu.log" 2>&1 \\
        || { cat "${sample_id}_emu.log" >&2; exit 1; }
    """

    stub:
    """
    cat > "${sample_id}_rel-abundance.tsv" <<'TSV'
tax_id	abundance	species	genus	family	order	class	phylum	clade	superkingdom	estimated counts
562	0.6	Escherichia coli	Escherichia	Enterobacteriaceae	Enterobacterales	Gammaproteobacteria	Pseudomonadota		Bacteria	540.0
1280	0.4	Staphylococcus aureus	Staphylococcus	Staphylococcaceae	Bacillales	Bacilli	Bacillota		Bacteria	360.0
unmapped	0								50
mapped_filtered	0								20
mapped_unclassified	0								30
TSV
    touch "${sample_id}_emu.log"
    """
}
