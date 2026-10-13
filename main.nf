#!/usr/bin/env nextflow
/*
========================================================================================
    16StJo - PIPELINE NANOPORE 16S (EMU) v1.0.0
    Compatible EPI2ME | Nextflow DSL2 | Docker
========================================================================================
    WORKFLOW :
       FASTQ (fichiers ou dossiers barcodeXX)
         → NANOPLOT (reads bruts)
         → CHOPPER (filtrage qualité / longueur)
             → NANOPLOT (reads filtrés)
             → EMU_ABUNDANCE (× chaque base de données choisie)
                 → EMU_COMBINE (tables multi-échantillons par base)
                 → EMU_REPORT  (rapport HTML interactif, comparaison des bases)
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMÈTRES DU PIPELINE
// ─────────────────────────────────────────────────────────────────────────────
params.fastq_dir      = null        // Dossier FASTQ (obligatoire) : fichiers ou sous-dossiers barcodeXX
params.sample_sheet   = null        // Optionnel : CSV "barcode,alias" pour renommer les barcodes
params.include_unclassified = false // Inclure le sous-dossier "unclassified" de MinKNOW
params.outdir         = "output"    // Dossier de sortie des rapports HTML
params.resultsdir     = "results"   // Dossier de sortie des fichiers intermédiaires

// Filtrage chopper (valeurs adaptées au 16S complet, amplicon ~1 500 pb)
params.min_quality    = 10          // Qualité moyenne minimale des reads (Q-score)
params.min_length     = 1200        // Longueur minimale des reads (pb)
params.max_length     = 1800        // Longueur maximale des reads (pb)
params.chopper_args   = ""          // Options supplémentaires passées à chopper
params.min_reads      = 50          // Nb minimal de reads filtrés pour lancer Emu sur un échantillon

// Emu
params.emu_db_dir     = "/data/emu_db"  // Dossier contenant les bases Emu (une par sous-dossier)
params.emu_dbs        = "emu"           // Base(s) à utiliser, séparées par des virgules (ex: "emu,silva,rdp")
params.emu_min_abundance = 0.0001       // Abondance minimale d'une espèce (--min-abundance)
params.emu_N          = 50              // Nb max d'alignements secondaires minimap2 (--N)
params.emu_args       = ""              // Options supplémentaires passées à "emu abundance"
params.emu_combine_ranks = "species,genus" // Rangs des tables combinées multi-échantillons

// Rapport
params.report_title   = "Analyse 16S Nanopore - Emu"
params.report_top_n   = 15          // Nb de taxons affichés dans les graphiques de composition

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// Chaque module = un outil = un fichier .nf dans le dossier modules/.
// ─────────────────────────────────────────────────────────────────────────────
include { NANOPLOT as NANOPLOT_RAW      } from './modules/nanoplot.nf'      // QC des reads bruts
include { NANOPLOT as NANOPLOT_FILTERED } from './modules/nanoplot.nf'      // QC des reads filtrés
include { CHOPPER       } from './modules/chopper.nf'       // Filtrage qualité/longueur
include { EMU_ABUNDANCE } from './modules/emu_abundance.nf' // Classification taxonomique
include { EMU_COMBINE   } from './modules/emu_combine.nf'   // Tables multi-échantillons
include { EMU_REPORT    } from './modules/emu_report.nf'    // Rapport HTML interactif

// ─────────────────────────────────────────────────────────────────────────────
// BANNIÈRE DE DÉMARRAGE
// Affiche les paramètres utilisés pour traçabilité dans les logs.
// ─────────────────────────────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════╗
║         16StJo - PIPELINE 16S NANOPORE (EMU)  v1.0.0     ║
╚══════════════════════════════════════════════════════════╝
  Dossier FASTQ   : ${params.fastq_dir}
  Sample sheet    : ${params.sample_sheet ?: 'aucune'}
  Dossier sortie  : ${params.outdir}
  Qualité min.    : ${params.min_quality}
  Longueur        : ${params.min_length} - ${params.max_length} bp
  Reads min.      : ${params.min_reads}
  Dossier bases   : ${params.emu_db_dir}
  Bases Emu       : ${params.emu_dbs}
  Abondance min.  : ${params.emu_min_abundance}
──────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
// VALIDATION DES PARAMÈTRES OBLIGATOIRES
// On arrête le pipeline proprement si un paramètre requis est absent.
// ─────────────────────────────────────────────────────────────────────────────
if (!params.fastq_dir) {
    error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
}
if (!params.emu_dbs) {
    error "ERREUR : --emu_dbs est obligatoire (ex : --emu_dbs emu,silva)."
}
if (params.min_length >= params.max_length) {
    error "ERREUR : --min_length (${params.min_length}) doit être inférieur à --max_length (${params.max_length})."
}

// Rangs acceptés par "emu combine-outputs"
def EMU_RANKS = ['tax_id', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
def combine_ranks = params.emu_combine_ranks.toString().tokenize(',')*.trim().findAll { it }
combine_ranks.each { rank ->
    if (!(rank in EMU_RANKS)) {
        error "ERREUR : rang '${rank}' invalide dans --emu_combine_ranks. Valeurs possibles : ${EMU_RANKS.join(', ')}"
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FONCTIONS UTILITAIRES
// ─────────────────────────────────────────────────────────────────────────────
def isFastq(f) {
    def name = f.name
    return f.isFile() && ['.fastq', '.fastq.gz', '.fq', '.fq.gz'].any { name.endsWith(it) }
}

def fastqName(f) {
    return f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
}

// Liste les échantillons du dossier d'entrée.
// LEÇON : deux organisations sont acceptées, comme dans EPI2ME :
//   1. fastq_dir/barcode01/*.fastq.gz  → un échantillon par sous-dossier
//      (sortie "fastq_pass" de MinKNOW, tous les fichiers sont concaténés)
//   2. fastq_dir/echantillon.fastq.gz  → un échantillon par fichier
def listSamples(fastq_dir, include_unclassified) {
    def dir = file(fastq_dir)
    if (!dir.exists() || !dir.isDirectory()) {
        error "ERREUR : --fastq_dir '${fastq_dir}' n'existe pas ou n'est pas un dossier."
    }
    def subdirs = dir.listFiles()
        .findAll { it.isDirectory() && it.listFiles().any { f -> isFastq(f) } }
        .findAll { include_unclassified || it.name != 'unclassified' }
        .sort { it.name }
    if (subdirs) {
        return subdirs.collect { d -> [d.name, d.listFiles().findAll { f -> isFastq(f) }.sort { it.name }] }
    }
    def files = dir.listFiles().findAll { isFastq(it) }.sort { it.name }
    if (!files) {
        error "ERREUR : aucun fichier FASTQ (.fastq, .fq, .fastq.gz, .fq.gz) trouvé dans '${fastq_dir}'."
    }
    return files.collect { f -> [fastqName(f), [f]] }
}

// Lit la sample sheet (CSV avec en-tête, colonnes "barcode" et "alias").
// Seuls les barcodes listés sont analysés ; ils sont renommés avec leur alias.
def readSampleSheet(path) {
    def lines = file(path, checkIfExists: true).readLines().findAll { it.trim() }
    def header = lines[0].split(',')*.trim()*.toLowerCase()
    def i_barcode = header.indexOf('barcode')
    def i_alias   = header.indexOf('alias')
    if (i_barcode < 0 || i_alias < 0) {
        error "ERREUR : la sample sheet doit contenir les colonnes 'barcode' et 'alias' (en-tête : ${lines[0]})."
    }
    def mapping = [:]
    lines.drop(1).each { line ->
        def cols = line.split(',', -1)*.trim()
        mapping[cols[i_barcode]] = cols[i_alias]
    }
    if (mapping.values().toSet().size() != mapping.size()) {
        error "ERREUR : la sample sheet contient des alias en double."
    }
    return mapping
}

// Résout les bases Emu : un nom simple = sous-dossier de --emu_db_dir,
// un chemin contenant "/" = dossier utilisé tel quel.
// Chaque base doit contenir species_taxid.fasta et taxonomy.tsv.
def resolveEmuDbs(db_list, db_dir) {
    def dbs = db_list.toString().tokenize(',')*.trim().findAll { it }.collect { entry ->
        def db = entry.contains('/') ? file(entry) : file("${db_dir}/${entry}")
        if (!db.isDirectory()) {
            error "ERREUR : base Emu '${entry}' introuvable (${db})."
        }
        ['species_taxid.fasta', 'taxonomy.tsv'].each { required ->
            if (!db.resolve(required).exists()) {
                error "ERREUR : la base Emu '${db}' ne contient pas ${required}."
            }
        }
        [db.name, db]
    }
    def names = dbs.collect { it[0] }
    if (names.toSet().size() != names.size()) {
        error "ERREUR : deux bases Emu portent le même nom de dossier (${names.join(', ')})."
    }
    return dbs
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Création du channel d'entrée ──────────────────────────────────────────
    // LEÇON : chaque élément est un tuple [sample_id, [fichiers FASTQ]].
    // Un échantillon peut avoir plusieurs fichiers (dossier barcodeXX de MinKNOW).
    def samples = listSamples(params.fastq_dir, params.include_unclassified)

    if (params.sample_sheet) {
        def sheet = readSampleSheet(params.sample_sheet)
        def missing = sheet.keySet() - samples.collect { it[0] }
        if (missing) {
            log.warn "Barcodes de la sample sheet absents de --fastq_dir : ${missing.join(', ')}"
        }
        samples = samples
            .findAll { sheet.containsKey(it[0]) }
            .collect { [sheet[it[0]], it[1]] }
        if (!samples) {
            error "ERREUR : aucun barcode de la sample sheet ne correspond aux dossiers de --fastq_dir."
        }
    }

    ch_raw = Channel.fromList(samples)

    // ── Bases de données Emu ──────────────────────────────────────────────────
    // Un tuple [nom_base, dossier] par base. Chaque échantillon sera classé
    // contre chacune d'elles pour pouvoir comparer les résultats.
    ch_dbs = Channel.fromList(resolveEmuDbs(params.emu_dbs, params.emu_db_dir))

    // ── QC des reads bruts ─────────────────────────────────────────────────────
    NANOPLOT_RAW(ch_raw.map { sample_id, reads -> tuple(sample_id, 'raw', reads) })

    // ── Filtrage qualité / longueur ───────────────────────────────────────────
    // chopper conserve les reads de qualité ≥ min_quality et de longueur
    // comprise entre min_length et max_length (taille de l'amplicon 16S).
    CHOPPER(ch_raw)

    // LEÇON : .branch() sépare un channel en plusieurs sous-channels.
    // Les échantillons avec trop peu de reads (ex : témoins négatifs) restent
    // visibles dans le QC du rapport, mais ne sont pas classés par Emu.
    ch_filtered = CHOPPER.out.reads
        .branch { sample_id, reads, n_reads ->
            pass: n_reads.toInteger() >= params.min_reads
            fail: true
        }

    ch_filtered.fail.subscribe { sample_id, reads, n_reads ->
        log.warn "${sample_id} : ${n_reads} reads après filtrage (< ${params.min_reads}) → exclu de la classification Emu"
    }

    ch_reads = ch_filtered.pass.map { sample_id, reads, n_reads -> tuple(sample_id, reads) }

    // ── QC des reads filtrés ───────────────────────────────────────────────────
    NANOPLOT_FILTERED(ch_reads.map { sample_id, reads -> tuple(sample_id, 'filtered', reads) })

    // ── Classification taxonomique Emu ────────────────────────────────────────
    // LEÇON : .combine() fait le produit cartésien échantillons × bases :
    // avec 10 échantillons et 3 bases, Emu est lancé 30 fois en parallèle.
    EMU_ABUNDANCE(ch_reads.combine(ch_dbs))

    // ── Tables combinées par base ─────────────────────────────────────────────
    // LEÇON : .groupTuple() regroupe toutes les tables d'une même base.
    ch_emu_by_db = EMU_ABUNDANCE.out.abundance
        .map { db_name, sample_id, tsv -> tuple(db_name, tsv) }
        .groupTuple()

    EMU_COMBINE(ch_emu_by_db, combine_ranks.join(' '))

    // ── Rapport HTML ──────────────────────────────────────────────────────────
    // On collecte toutes les sorties : stats NanoPlot (bruts + filtrés),
    // nombre de reads par échantillon après chopper, et un dossier par base Emu.
    ch_nanostats = NANOPLOT_RAW.out.stats
        .mix(NANOPLOT_FILTERED.out.stats)
        .collect()

    ch_read_counts = CHOPPER.out.reads
        .map { sample_id, reads, n_reads -> "${sample_id}\t${n_reads}" }
        .collectFile(name: 'read_counts.tsv', newLine: true, sort: true)

    ch_emu_dirs = EMU_COMBINE.out.db_dir
        .collect()
        .ifEmpty([])

    def run_info = [
        pipeline          : "16StJo v${workflow.manifest.version}",
        date              : new Date().format('yyyy-MM-dd HH:mm'),
        fastq_dir         : params.fastq_dir.toString(),
        sample_sheet      : params.sample_sheet ? params.sample_sheet.toString() : 'aucune',
        min_quality       : params.min_quality,
        min_length        : params.min_length,
        max_length        : params.max_length,
        chopper_args      : params.chopper_args ?: '-',
        min_reads         : params.min_reads,
        emu_dbs           : params.emu_dbs.toString(),
        emu_min_abundance : params.emu_min_abundance,
        emu_N             : params.emu_N,
        emu_args          : params.emu_args ?: '-',
        resultsdir        : params.resultsdir.toString()
    ]

    EMU_REPORT(
        ch_nanostats,
        ch_read_counts,
        ch_emu_dirs,
        groovy.json.JsonOutput.toJson(run_info)
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// RÉSUMÉ DE FIN DE PIPELINE
// S'exécute automatiquement à la fin, succès ou échec.
// ─────────────────────────────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
╔══════════════════════════════════════════════════════════╗
║                  PIPELINE TERMINÉ !                      ║
╚══════════════════════════════════════════════════════════╝
  Statut    : ${workflow.success ? '✅ Succès' : '❌ Échec'}
  Durée     : ${workflow.duration}
  Rapports  : ${params.outdir}/
  Résultats : ${params.resultsdir}/
──────────────────────────────────────────────────────────
""".stripIndent()
}
