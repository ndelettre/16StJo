# 16StJo

Pipeline d'analyse taxonomique **16S Oxford Nanopore** basé sur [Emu](https://github.com/treangenlab/emu) — Nextflow DSL2, Docker, compatible EPI2ME.

## Description

Le pipeline filtre les reads 16S, contrôle leur qualité, les classe avec Emu contre **une ou plusieurs bases de données** (pour pouvoir les comparer) et produit un rapport HTML interactif.

### Workflow

```
FASTQ (fichiers ou dossiers barcodeXX)
  → NANOPLOT (reads bruts)
  → CHOPPER (qualité + longueur)
      → NANOPLOT (reads filtrés)
      → EMU_ABUNDANCE (× chaque base choisie)
          → EMU_COMBINE (tables multi-échantillons par base)
          → EMU_REPORT  (rapport HTML : QC, composition, diversité, comparaison des bases)
```

## Prérequis

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Docker](https://www.docker.com/)
- Une ou plusieurs bases Emu, une par sous-dossier de `/data/emu_db` (chaque base contient `species_taxid.fasta` et `taxonomy.tsv`)

### Téléchargement des bases Emu

Les bases pré-construites sont sur [OSF](https://osf.io/56uf7/) :

| Base | Nom OSF | Contenu |
|------|---------|---------|
| Emu (défaut) | `emu` | rrnDB v5.6 + NCBI 16S RefSeq (bactéries et archées) |
| RDP | `rdp` | RDP v11.5, taxonomie NCBI |
| SILVA | `silva` | SILVA v138.1 (base DADA2 au niveau espèce) |
| SILVA 138.2 | `silva-138.2` | non validée par les auteurs d'Emu |

```bash
pip install osfclient
for db in emu silva rdp; do
    mkdir -p /data/emu_db/$db
    cd /data/emu_db/$db
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/$db.tar
    tar -xvf $db.tar && rm $db.tar
done
```

Vérifier ensuite que `species_taxid.fasta` et `taxonomy.tsv` sont **directement** dans chaque dossier (`/data/emu_db/silva/taxonomy.tsv`) ; si l'archive a créé un sous-dossier, remonter ses fichiers d'un niveau. Le pipeline vérifie leur présence au démarrage.

## Installation

### Via EPI2ME Desktop

1. Ouvrir EPI2ME Desktop
2. **Workflows** → **Import workflow** → URL GitHub `https://github.com/ndelettre/16StJo`
3. Ou **Import Local Workflow** → sélectionner ce dossier

Les paramètres apparaissent dans l'interface grâce à `nextflow_schema.json`.

### En ligne de commande

```bash
nextflow run ndelettre/16StJo \
    --fastq_dir /chemin/vers/fastq_pass \
    --sample_sheet samples.csv \
    --emu_db_dir /data/emu_db \
    --emu_dbs emu,silva,rdp
```

## Entrées

Deux organisations du dossier `fastq_dir` sont acceptées :

```
fastq_pass/                    ← sortie MinKNOW : un sous-dossier par barcode
├── barcode01/*.fastq.gz       ← tous les fichiers d'un barcode sont concaténés
├── barcode02/*.fastq.gz
└── unclassified/              ← ignoré (sauf --include_unclassified true)

fastq/                         ← ou un fichier par échantillon
├── echantillon1.fastq.gz
└── echantillon2.fastq
```

**Sample sheet** (optionnelle, format EPI2ME) : CSV avec les colonnes `barcode` et `alias`. Seuls les barcodes listés sont analysés, et ils sont renommés avec leur alias.

```csv
barcode,alias
barcode01,Patient_A
barcode02,Patient_B
barcode03,Temoin_negatif
```

## Paramètres

### Entrées / Sorties

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `fastq_dir` | Dossier FASTQ (fichiers ou sous-dossiers barcodeXX) | **obligatoire** |
| `sample_sheet` | CSV `barcode,alias` | — |
| `include_unclassified` | Analyser le dossier `unclassified` de MinKNOW | `false` |
| `outdir` | Dossier du rapport HTML | `output` |
| `resultsdir` | Dossier des fichiers intermédiaires | `results` |

### Filtrage (chopper)

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `min_quality` | Qualité moyenne minimale (Q-score) | `10` |
| `min_length` | Longueur minimale (pb) | `1200` |
| `max_length` | Longueur maximale (pb) | `1800` |
| `chopper_args` | Options supplémentaires de chopper | — |
| `min_reads` | Reads minimum après filtrage pour lancer Emu | `50` |

Les valeurs par défaut conviennent au 16S complet (~1 500 pb, kit SQK-16S114). Pour un amplicon plus court (ex : V3-V4), adapter `min_length` et `max_length`.

Les échantillons sous `min_reads` (ex : témoins négatifs) restent dans le QC du rapport mais ne sont pas classés.

### Emu

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `emu_db_dir` | Dossier contenant les bases (une par sous-dossier) | `/data/emu_db` |
| `emu_dbs` | Base(s) à utiliser, séparées par des virgules : noms de sous-dossiers ou chemins complets | `emu` |
| `emu_min_abundance` | Abondance minimale d'une espèce (`--min-abundance`) | `0.0001` |
| `emu_N` | Alignements secondaires minimap2 max par read (`--N`) ; le baisser réduit la RAM | `50` |
| `emu_args` | Options supplémentaires de `emu abundance` (ex : `--min-pid 90`) | — |
| `emu_combine_ranks` | Rangs des tables combinées | `species,genus` |

Chaque échantillon est classé contre **chaque** base de `emu_dbs`. Exemple pour comparer trois bases :

```bash
--emu_dbs emu,silva,rdp
```

### Rapport

| Paramètre | Description | Défaut |
|-----------|-------------|--------|
| `report_title` | Titre du rapport | `Analyse 16S Nanopore - Emu` |
| `report_top_n` | Nombre de taxons affichés dans les graphiques | `15` |

## Sorties

```
output/
├── 16StJo_report.html               ← Rapport interactif (autonome, lisible hors ligne)
├── pipeline_report.html             ← Métriques d'exécution Nextflow
├── pipeline_timeline.html           ← Chronologie d'exécution
└── pipeline_dag.html                ← Graphe du pipeline

results/
├── chopper/                         ← Reads filtrés + bilan du filtrage
├── nanoplot/<échantillon>_raw/      ← Rapport NanoPlot des reads bruts
├── nanoplot/<échantillon>_filtered/ ← Rapport NanoPlot des reads filtrés
└── emu/<base>/
    ├── <échantillon>_rel-abundance.tsv       ← Abondances + comptages estimés
    ├── emu-combined-<rang>.tsv               ← Abondances, tous échantillons
    └── emu-combined-<rang>-counts.tsv        ← Comptages, tous échantillons
```

### Rapport HTML (version 1)

- **Résumé** : nombre d'échantillons, de reads, bases utilisées, paramètres du run
- **Contrôle qualité** : reads avant / après chopper, longueur et qualité médianes, N50, échantillons exclus
- **Classification Emu** : part des reads classés / non alignés / filtrés, par base
- **Composition taxonomique** : barres empilées, base et rang (espèce → phylum) au choix
- **Diversité alpha** : richesse, Shannon, Simpson inverse
- **Comparaison des bases** : heatmap taxons × (échantillon | base), dissimilarité de Bray-Curtis entre bases, taxon dominant par base

Le contenu du rapport est généré par `bin/emu_report.py` ; chaque section est une fonction `section_*` facile à modifier ou à compléter. Pour un PDF : ouvrir le rapport dans un navigateur → Imprimer → Enregistrer en PDF (toutes les vues sont alors affichées).

## Logiciels utilisés

| Outil | Version | Conteneur |
|-------|---------|-----------|
| chopper | 0.12.0b | `community.wave.seqera.io/library/chopper:0.12.0b--f39108dd84394289` |
| NanoPlot | 1.47.0 | `quay.io/biocontainers/nanoplot:1.47.0--pyhdfd78af_0` |
| Emu | 3.6.2 | `community.wave.seqera.io/library/emu:3.6.2--0adcfa4809f061f6` |
| Rapport (pandas, plotly) | — | `nicolasdelettre/phyloreport:1.0.1` |

Les conteneurs chopper, NanoPlot et Emu sont ceux des modules nf-core.

## Notes

- Les noms de taxons diffèrent parfois d'une base à l'autre (ex : *Escherichia-Shigella* dans SILVA, *Escherichia* dans NCBI). Le rapport compare les noms tels quels : ces différences de nomenclature apparaissent comme des différences de composition.
- La RAM d'Emu dépend surtout de la taille de la base. En cas d'erreur mémoire, diminuer `emu_N` ou augmenter `memory` du label `process_medium` dans `nextflow.config`.
- Test rapide de la logique du pipeline sans Docker ni données réelles : `nextflow run main.nf -stub -profile local --fastq_dir <dossier> --emu_db_dir <dossier>`.
- En cas d'erreur, consulter `work/xx/yyyy/.command.err` (chemin indiqué par Nextflow).
