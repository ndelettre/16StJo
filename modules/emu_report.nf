/*
========================================================================================
    MODULE : EMU_REPORT — Rapport HTML interactif
========================================================================================
    Rassemble le QC (NanoPlot + chopper) et les résultats Emu de toutes les
    bases dans un seul fichier HTML autonome (consultable hors ligne) :
      - QC des reads avant / après filtrage
      - taux de classification par base
      - composition taxonomique (barres empilées, rang au choix)
      - diversité alpha
      - comparaison des bases entre elles (si plusieurs bases)

    Script : bin/emu_report.py (pandas + plotly)
========================================================================================
*/

process EMU_REPORT {

    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path nanostats,   stageAs: 'nanostats/*'
    path read_counts
    path emu_dirs,    stageAs: 'emu/*'
    val  run_info

    output:
    path "16StJo_report.html", emit: html

    script:
    // LEÇON : le JSON des paramètres est écrit via un heredoc entre quotes
    // ('JSON') pour que bash n'interprète ni les $ ni les guillemets.
    """
    set -euo pipefail

    cat > run_info.json <<'JSON'
${run_info}
JSON

    python3 ${projectDir}/bin/emu_report.py \\
        --nanostats_dir nanostats \\
        --read_counts   ${read_counts} \\
        --emu_dir       emu \\
        --run_info      run_info.json \\
        --title         "${params.report_title}" \\
        --top_n         ${params.report_top_n} \\
        --min_reads     ${params.min_reads} \\
        --output        16StJo_report.html
    """

    stub:
    """
    touch 16StJo_report.html
    """
}
