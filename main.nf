#!/usr/bin/env nextflow
/*
========================================================================================
    ARN16StJo — PIPELINE 16S FULL-LENGTH NANOPORE (EMU)  v0.2.0
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW :
       NANOSTAT (reads bruts)
       CHOPPER → NANOSTAT (reads filtrés)
               → HOST_MAP + HOST_FILTER (comptage et retrait des reads de l'hôte)
               → EMU_ABUNDANCE → EMU_COMBINE (tables espèce / genre / famille)
                               → SPIKE_IN (reads du contrôle interne Zymo)

    Étapes à venir : rapport HTML, mode recherche (diversité).
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
include { HOST_MAP                      } from './modules/host_map.nf'      // Alignement sur l'hôte
include { HOST_FILTER                   } from './modules/host_filter.nf'   // Comptage/retrait hôte
include { EMU_ABUNDANCE                 } from './modules/emu_abundance.nf' // Classification EMU
include { EMU_COMBINE                   } from './modules/emu_combine.nf'   // Tables multi-échantillons
include { SPIKE_IN                      } from './modules/spike_in.nf'      // Contrôle interne Zymo

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // LEÇON : la bannière et les validations sont DANS le bloc workflow :
    // la syntaxe stricte de Nextflow (par défaut dans les futures versions)
    // interdit les instructions en dehors d'un workflow ou d'un process.

    // ── Génome hôte (optionnel) ───────────────────────────────────────────────
    // LEÇON : un paramètre laissé vide arrive comme "" (EPI2ME) ou comme true
    // (ligne de commande "--host_reference" sans valeur). On ne garde le chemin
    // que s'il s'agit d'un texte non vide ; sinon l'étape hôte est ignorée.
    def host_reference = (params.host_reference instanceof CharSequence && params.host_reference.toString().trim())
        ? params.host_reference.toString().trim()
        : null

    // ── Bannière de démarrage ─────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║     ARN16StJo - 16S FULL-LENGTH NANOPORE (EMU)  v0.2.0   ║
    ╚══════════════════════════════════════════════════════════╝
      Dossier FASTQ   : ${params.fastq_dir}
      Dossier sortie  : ${params.outdir}
      Longueur reads  : ${params.min_length} - ${params.max_length} bp
      Qualité min.    : Q${params.min_quality}
      Génome hôte     : ${host_reference ?: 'aucun (étape ignorée)'}
      Base EMU        : ${params.emu_db}
      Contrôle interne: ${params.spike_in}
    ──────────────────────────────────────────────────────────
    """.stripIndent()

    // ── Validation des paramètres ─────────────────────────────────────────────
    if (!params.fastq_dir) {
        error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
    }
    if (params.min_length >= params.max_length) {
        error "ERREUR : --min_length (${params.min_length}) doit être inférieur à --max_length (${params.max_length})."
    }
    if (host_reference && !(host_reference =~ /\.((fa|fasta|fna)(\.gz)?|mmi)$/)) {
        error "ERREUR : --host_reference doit être un FASTA (.fa, .fasta, .fna, éventuellement .gz) ou un index minimap2 (.mmi)."
    }
    if (host_reference && !file(host_reference).exists()) {
        error "ERREUR : génome hôte introuvable : ${host_reference}\n(voir le README pour le télécharger, ou laisser --host_reference vide pour ignorer l'étape)."
    }
    ['species_taxid.fasta', 'taxonomy.tsv'].each { db_file ->
        if (!file("${params.emu_db}/${db_file}").exists()) {
            error "ERREUR : ${db_file} introuvable dans la base EMU --emu_db (${params.emu_db}).\n(voir le README pour télécharger la base)."
        }
    }

    // ── Contrôle interne ZymoBIOMICS ──────────────────────────────────────────
    // Souches attendues pour chaque contrôle (fiches Zymo D6320 et D6321).
    def spike_in_controls = [
        spike_in_I : [name   : 'Spike-in Control I (D6320, High Microbial Load)',
                      species: ['Imtechella halotolerans', 'Allobacillus halotolerans']],
        spike_in_II: [name   : 'Spike-in Control II (D6321, Low Microbial Load)',
                      species: ['Truepera radiovictrix', 'Imtechella halotolerans', 'Allobacillus halotolerans']]
    ]
    if (params.spike_in != 'aucun' && !spike_in_controls.containsKey(params.spike_in)) {
        error "ERREUR : --spike_in doit valoir 'aucun', 'spike_in_I' ou 'spike_in_II' (reçu : '${params.spike_in}')."
    }
    def spike_in = spike_in_controls[params.spike_in]

    // Vérifie que les souches du contrôle existent dans la base EMU : sinon,
    // leurs reads seraient attribués à une espèce proche et fausseraient tout.
    if (spike_in) {
        def db_species = file("${params.emu_db}/taxonomy.tsv").readLines().collect { line -> line.tokenize('\t')[1] } as Set
        spike_in.species.findAll { sp -> !db_species.contains(sp) }.each { sp ->
            log.warn "Contrôle interne : '${sp}' est absente de la base EMU ${params.emu_db}"
        }
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

    // ── Reads de l'hôte : comptage et retrait ─────────────────────────────────
    // HOST_MAP aligne les reads de TOUS les échantillons en un seul passage
    // (le génome humain n'est indexé/chargé qu'une fois), puis HOST_FILTER
    // retire les reads hôte de chaque échantillon.
    //
    // LEÇON : .collect() rassemble les FASTQ de tous les échantillons en une
    // seule liste. Comme toutes les entrées de HOST_MAP sont alors des "value
    // channels", sa sortie host_reads.tsv en est une aussi : elle peut être
    // réutilisée par HOST_FILTER pour chacun des échantillons.
    if (host_reference) {
        ch_host_reference = channel.fromPath(host_reference, checkIfExists: true).first()
        HOST_MAP(CHOPPER.out.reads.map { _sample_id, reads -> reads }.collect(), ch_host_reference)
        HOST_FILTER(CHOPPER.out.reads, HOST_MAP.out.hits)
        ch_reads_emu = HOST_FILTER.out.reads

        // Table récapitulative : une ligne par échantillon
        HOST_FILTER.out.stats
            .map { _sample_id, stats -> stats }
            .collectFile(name: 'host_summary.tsv', keepHeader: true, skip: 1,
                         sort: { stats -> stats.name }, storeDir: "${params.resultsdir}/host")
    } else {
        log.warn "Aucun génome hôte (--host_reference) : les reads de l'hôte ne sont ni comptés ni retirés."
        ch_reads_emu = CHOPPER.out.reads
    }

    // ── Classification taxonomique EMU ────────────────────────────────────────
    ch_emu_db = channel.fromPath(params.emu_db, type: 'dir', checkIfExists: true).first()
    EMU_ABUNDANCE(ch_reads_emu, ch_emu_db)

    ch_abundance_tables = EMU_ABUNDANCE.out.abundance
        .map { _sample_id, table -> table }
        .collect()

    // Tables multi-échantillons (espèce, genre, famille)
    EMU_COMBINE(ch_abundance_tables)

    // ── Contrôle interne ──────────────────────────────────────────────────────
    if (spike_in) {
        SPIKE_IN(ch_abundance_tables, spike_in.name, spike_in.species.join(','))
    }

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
