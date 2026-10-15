# ARN16StJo

Pipeline Nextflow d'analyse du **16S full-length** séquencé sur **Oxford Nanopore** (R10.4.1), avec classification taxonomique par [EMU](https://github.com/treangenlab/emu) — compatible **EPI2ME Desktop**.

> 🚧 Pipeline en cours de développement. Version actuelle : **0.1.0** (étape 1 : contrôle qualité).

## Description

Le pipeline prend en entrée des reads déjà basecallés et démultiplexés (un FASTQ par échantillon, amorces et adaptateurs retirés), par exemple issus du *Microbial Amplicon Barcoding Kit 24 V14*.

### Workflow

```
FASTQ ─┬─ NANOSTAT (reads bruts)
       └─ CHOPPER (longueur + qualité) ── NANOSTAT (reads filtrés)
```

Étapes prévues :

| Étape | Contenu | Statut |
|---|---|---|
| 1 | Squelette + QC : NanoStat → Chopper → NanoStat | ✅ |
| 2 | Comptage et retrait des reads de l'hôte (minimap2, GRCh38 par défaut) | ⏳ |
| 3 | Classification EMU (base rrnDB 5.10 + NCBI mars 2026 par défaut) | ⏳ |
| 4 | Rapport HTML unique (QC, hôte, contrôle interne Zymo, abondances) | ⏳ |
| 5 | Mode recherche optionnel (diversité alpha/bêta, export phyloseq) | ⏳ |

## Prérequis

- [Nextflow](https://www.nextflow.io/) ≥ 22.10.0
- [Docker](https://www.docker.com/)

## Installation

### Via EPI2ME Desktop

1. Ouvrir EPI2ME Desktop
2. **Workflows** → **Import workflow** → `https://github.com/ndelettre/arn16stjo`
3. Ou **Import Local Workflow** → sélectionner ce dossier

> ⚠️ Le nom du dépôt et `manifest.name` doivent rester **en minuscules** : EPI2ME Desktop
> cherche le workflow dans un dossier au nom converti en minuscules, alors que `nextflow pull`
> garde la casse de l'URL. Avec des majuscules, l'import échoue (« Error refreshing workflow config »).

### En ligne de commande

```bash
nextflow run ndelettre/arn16stjo \
    --fastq_dir /chemin/vers/fastq
```

## Données d'entrée

Un dossier contenant **un fichier FASTQ par échantillon** (`.fastq`, `.fq`, `.fastq.gz` ou `.fq.gz`).
Le nom de l'échantillon est le nom du fichier sans l'extension :

| Fichier | Échantillon |
|---|---|
| `patient01.fastq.gz` | `patient01` |
| `temoin_neg.fq` | `temoin_neg` |
| `run2.barcode05.fastq` | `run2.barcode05` |

Les témoins (négatif, positif) sont traités comme des échantillons.
Deux fichiers donnant le même nom d'échantillon (ex. `S1.fastq` et `S1.fastq.gz`) arrêtent le pipeline avec un message d'erreur.

## Paramètres

| Paramètre | Défaut | Description |
|---|---|---|
| `--fastq_dir` | *(obligatoire)* | Dossier contenant les FASTQ |
| `--min_length` | `1200` | Longueur minimale des reads (bp) |
| `--max_length` | `1800` | Longueur maximale des reads (bp) |
| `--min_quality` | `10` | Q-score moyen minimal des reads |
| `--outdir` | `output` | Dossier des rapports |
| `--resultsdir` | `results` | Dossier des résultats intermédiaires |

## Sorties

```
results/
├── nanostat/
│   ├── {échantillon}_raw_nanostat.tsv        # QC avant filtrage
│   └── {échantillon}_filtered_nanostat.tsv   # QC après filtrage
└── chopper/
    ├── {échantillon}_filtered.fastq.gz       # reads conservés
    └── {échantillon}_chopper.log             # "Kept X reads out of Y reads"
output/
├── pipeline_report.html                      # ressources utilisées par tâche
└── pipeline_timeline.html                    # chronologie d'exécution
```

Un échantillon sans read (barcode vide, ou tous les reads écartés par Chopper) n'arrête pas le run : un avertissement est affiché dans le log et ses statistiques NanoStat sont à 0.

## Images Docker

Toutes les images sont figées sur une version exacte pour garantir la reproductibilité.

| Outil | Image |
|---|---|
| NanoStat | `quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0` |
| Chopper | `quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0` |
