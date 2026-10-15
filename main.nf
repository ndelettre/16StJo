#!/usr/bin/env nextflow
/*
========================================================================================
    ARN16StJo — PIPELINE 16S FULL-LENGTH NANOPORE (EMU)  v0.1.0
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW (étape 1 : contrôle qualité) :
       NANOSTAT (reads bruts) → CHOPPER → NANOSTAT (reads filtrés)

    Étapes à venir : retrait des reads de l'hôte, classification EMU,
    rapport HTML, mode recherche (diversité).
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// LEÇON : un même process peut être importé deux fois sous deux noms
// différents avec "as". Chaque alias est un process indépendant : on peut
// donc lancer NanoStat avant ET après Chopper sans dupliquer le module.
// ─────────────────────────────────────────────────────────────────────────────
include { NANOSTAT as NANOSTAT_RAW      } from './modules/nanostat.nf' // QC reads bruts
include { CHOPPER                       } from './modules/chopper.nf'  // Filtrage longueur/qualité
include { NANOSTAT as NANOSTAT_FILTERED } from './modules/nanostat.nf' // QC reads filtrés

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // LEÇON : la bannière et les validations sont DANS le bloc workflow :
    // la syntaxe stricte de Nextflow (par défaut dans les futures versions)
    // interdit les instructions en dehors d'un workflow ou d'un process.

    // ── Bannière de démarrage ─────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║     ARN16StJo - 16S FULL-LENGTH NANOPORE (EMU)  v0.1.0   ║
    ╚══════════════════════════════════════════════════════════╝
      Dossier FASTQ   : ${params.fastq_dir}
      Dossier sortie  : ${params.outdir}
      Longueur reads  : ${params.min_length} - ${params.max_length} bp
      Qualité min.    : Q${params.min_quality}
    ──────────────────────────────────────────────────────────
    """.stripIndent()

    // ── Validation des paramètres ─────────────────────────────────────────────
    if (!params.fastq_dir) {
        error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
    }
    if (params.min_length >= params.max_length) {
        error "ERREUR : --min_length (${params.min_length}) doit être inférieur à --max_length (${params.max_length})."
    }

    // ── Création du channel d'entrée ──────────────────────────────────────────
    // Un fichier FASTQ = un échantillon. Le nom de l'échantillon est le nom du
    // fichier sans l'extension .fastq/.fq(.gz) : "patient.01.fastq.gz" → "patient.01".
    // LEÇON : on n'utilise pas getSimpleName() car il coupe au PREMIER point
    // ("patient.01.fastq.gz" deviendrait "patient").
    //
    // LEÇON : groupTuple() regroupe les fichiers par nom d'échantillon. Si deux
    // fichiers donnent le même nom (ex : S1.fastq et S1.fastq.gz), on arrête
    // le pipeline au lieu d'écraser silencieusement un des deux.
    ch_fastq = channel
        .fromPath(
            "${params.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}",
            checkIfExists: true
        )
        .map { file ->
            def sample_id = file.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample_id, file)
        }
        .groupTuple()
        .map { sample_id, files ->
            if (files.size() > 1) {
                error "ERREUR : plusieurs fichiers pour l'échantillon '${sample_id}' : ${files*.name.join(', ')}"
            }
            tuple(sample_id, files[0])
        }

    // ── Contrôle qualité avant filtrage ───────────────────────────────────────
    // LEÇON : le deuxième argument ("raw") est une simple valeur passée au
    // process : elle sert à nommer les fichiers et le dossier de sortie.
    NANOSTAT_RAW(ch_fastq, "raw")

    // ── Filtrage longueur / qualité ───────────────────────────────────────────
    // Chopper garde les reads entre min_length et max_length avec un Q-score
    // moyen ≥ min_quality (16S full-length ≈ 1 500 bp).
    CHOPPER(ch_fastq)

    // Alerte dans le log EPI2ME pour chaque échantillon sans read après filtrage
    // (le run continue : l'échantillon apparaîtra avec 0 read dans le rapport).
    CHOPPER.out.log
        .filter { _sample_id, chopper_log -> chopper_log.text =~ /^Kept 0 reads/ }
        .subscribe { sample_id, chopper_log ->
            log.warn "${sample_id} : aucun read conservé par Chopper (${chopper_log.text.trim()})"
        }

    // ── Contrôle qualité après filtrage ───────────────────────────────────────
    NANOSTAT_FILTERED(CHOPPER.out.reads, "filtered")

    // ── Résumé de fin de pipeline (succès ou échec) ───────────────────────────
    // LEÇON : à l'intérieur du bloc onComplete, les noms `workflow` et `params`
    // ne sont plus accessibles : on en garde des copies locales avant.
    def run        = workflow
    def outdir     = params.outdir
    def resultsdir = params.resultsdir
    run.onComplete {
        log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║                  PIPELINE TERMINÉ !                      ║
    ╚══════════════════════════════════════════════════════════╝
      Statut    : ${run.success ? '✅ Succès' : '❌ Échec'}
      Durée     : ${run.duration}
      Rapports  : ${outdir}/
      Résultats : ${resultsdir}/
    ──────────────────────────────────────────────────────────
    """.stripIndent()
    }
}
