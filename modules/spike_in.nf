/*
========================================================================================
    MODULE : SPIKE_IN — Reads du contrôle interne ZymoBIOMICS
========================================================================================
    Compte les reads attribués par EMU aux souches du contrôle interne, pour
    chaque échantillon (script bin/spike_in_counts.py).

    Souches attendues (définies dans main.nf) :
    - Spike-in Control I  (D6320) : Imtechella halotolerans, Allobacillus halotolerans
    - Spike-in Control II (D6321) : Truepera radiovictrix, Imtechella halotolerans,
                                    Allobacillus halotolerans

    Sorties :
    - spike_in_counts.tsv → sample | control | species | expected |
                            estimated_reads | pct_classified | pct_total
========================================================================================
*/

process SPIKE_IN {

    label 'process_mono'

    publishDir "${params.resultsdir}/emu", mode: 'copy'

    input:
    path abundance_tables   // tous les {sample_id}_rel-abundance.tsv
    val control_name        // ex. "Spike-in Control I (D6320)"
    val species             // espèces attendues, séparées par des virgules

    output:
    path "spike_in_counts.tsv", emit: counts

    script:
    """
    spike_in_counts.py \\
        --control "${control_name}" \\
        --species "${species}" \\
        --output spike_in_counts.tsv \\
        *_rel-abundance.tsv
    """
}
