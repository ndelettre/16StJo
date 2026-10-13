#!/usr/bin/env python3
"""
========================================================================================
    emu_report.py — Rapport HTML interactif du pipeline 16StJo
========================================================================================
Rassemble dans un seul fichier HTML autonome (plotly.js intégré, lisible hors
ligne) :
  1. le QC des reads (NanoPlot avant / après chopper)
  2. le taux de classification Emu par base de données
  3. la composition taxonomique (barres empilées, base et rang au choix)
  4. la diversité alpha (richesse, Shannon, Simpson inverse)
  5. la comparaison des bases entre elles (heatmap, Bray-Curtis, taxon dominant)

Entrées :
  --nanostats_dir : fichiers {sample}_{raw|filtered}_NanoStats.txt
  --read_counts   : TSV "sample<TAB>nb reads après chopper"
  --emu_dir       : un sous-dossier par base, contenant les {sample}_rel-abundance.tsv
  --run_info      : JSON des paramètres du run

Dépendances : pandas, numpy, plotly
========================================================================================
"""

import argparse
import html
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

# Lignes spéciales ajoutées par Emu à la fin de chaque table
EMU_SPECIAL_ROWS = ["unmapped", "mapped_filtered", "mapped_unclassified"]

# Rangs proposés dans le rapport, du plus fin au plus large
RANKS = ["species", "genus", "family", "order", "class", "phylum"]
RANK_LABELS = {
    "species": "Espèce", "genus": "Genre", "family": "Famille",
    "order": "Ordre", "class": "Classe", "phylum": "Phylum",
}

OTHER_LABEL = "Autres"
UNKNOWN_LABEL = "Non assigné"

# Palette qualitative (Plotly + D3 + Pastel) pour colorer les taxons
PALETTE = [
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692",
    "#B6E880", "#FF97FF", "#FECB52", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
    "#9467BD", "#8C564B", "#E377C2", "#BCBD22", "#17BECF", "#66C5CC", "#F6CF71",
    "#F89C74", "#DCB0F2", "#87C55F", "#9EB9F3", "#FE88B1", "#C9DB74", "#8BE0A4",
    "#B497E7", "#D3B484",
]
OTHER_COLOR = "#BBBBBB"
UNKNOWN_COLOR = "#777777"

# Couleurs des catégories de classification Emu
CLASSIF_COLORS = {
    "Classés": "#2CA02C",
    "Alignés non classés": "#FFA15A",
    "Alignés filtrés": "#FECB52",
    "Non alignés": "#D62728",
}


# ─────────────────────────────────────────────────────────────────────────────
# OUTILS
# ─────────────────────────────────────────────────────────────────────────────
def natural_key(text):
    """Tri naturel : barcode2 avant barcode10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text))]


def esc(value):
    return html.escape(str(value))


def fmt_int(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "-"
    return f"{int(round(value)):,}".replace(",", " ")


def fmt_float(value, digits=1):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "-"
    return f"{value:.{digits}f}"


def fig_html(fig, div_id):
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id=div_id,
                       config={"displaylogo": False, "responsive": True})


def slug(text):
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(text))


# ─────────────────────────────────────────────────────────────────────────────
# LECTURE DES DONNÉES
# ─────────────────────────────────────────────────────────────────────────────
def parse_nanostats(path):
    """Lit un NanoStats.txt (format TSV --tsv_stats ou format texte classique)."""
    stats = {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if "\t" in line:
                key, _, value = line.partition("\t")
            elif ":" in line:
                key, _, value = line.partition(":")
            else:
                continue
            key = key.strip().lower().replace(" ", "_").replace(",", "")
            value = value.strip().split()[0].replace(",", "") if value.strip() else ""
            try:
                stats[key] = float(value)
            except ValueError:
                continue
    # Harmonise quelques noms entre versions de NanoPlot
    aliases = {
        "number_of_reads": ["number_of_reads"],
        "number_of_bases": ["number_of_bases", "total_bases"],
        "median_read_length": ["median_read_length"],
        "mean_read_length": ["mean_read_length"],
        "n50": ["n50", "read_length_n50"],
        "median_qual": ["median_qual", "median_read_quality"],
        "mean_qual": ["mean_qual", "mean_read_quality"],
    }
    return {name: next((stats[k] for k in keys if k in stats), np.nan)
            for name, keys in aliases.items()}


def load_qc(nanostats_dir, read_counts_path):
    """Table QC : une ligne par échantillon, stats brutes et filtrées."""
    rows = {}
    for path in sorted(Path(nanostats_dir).glob("*_NanoStats.txt")):
        match = re.match(r"(.+)_(raw|filtered)_NanoStats\.txt$", path.name)
        if not match:
            continue
        sample, stage = match.groups()
        for key, value in parse_nanostats(path).items():
            rows.setdefault(sample, {})[f"{stage}_{key}"] = value

    counts = {}
    if read_counts_path and os.path.exists(read_counts_path):
        with open(read_counts_path) as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) == 2 and parts[1].strip().isdigit():
                    counts[parts[0]] = int(parts[1])

    for sample in set(rows) | set(counts):
        rows.setdefault(sample, {})["kept_reads"] = counts.get(sample, np.nan)

    qc = pd.DataFrame.from_dict(rows, orient="index")
    if qc.empty:
        return qc
    qc = qc.loc[sorted(qc.index, key=natural_key)]
    if "raw_number_of_reads" in qc and "kept_reads" in qc:
        qc["pct_kept"] = 100 * qc["kept_reads"] / qc["raw_number_of_reads"]
    return qc


def read_emu_table(path):
    df = pd.read_csv(path, sep="\t", dtype={"tax_id": str})
    df["tax_id"] = df["tax_id"].astype(str)
    df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce").fillna(0.0)
    if "estimated counts" in df.columns:
        df["estimated counts"] = pd.to_numeric(df["estimated counts"], errors="coerce").fillna(0.0)
    return df


def load_emu(emu_dir):
    """Retourne {base: {échantillon: DataFrame Emu}}."""
    data = {}
    emu_dir = Path(emu_dir)
    if not emu_dir.is_dir():
        return data
    for db_dir in sorted(p for p in emu_dir.iterdir() if p.is_dir()):
        tables = {}
        for path in db_dir.glob("*_rel-abundance.tsv"):
            sample = path.name[: -len("_rel-abundance.tsv")]
            tables[sample] = read_emu_table(path)
        if tables:
            data[db_dir.name] = dict(sorted(tables.items(), key=lambda kv: natural_key(kv[0])))
    return data


# ─────────────────────────────────────────────────────────────────────────────
# CALCULS
# ─────────────────────────────────────────────────────────────────────────────
def taxa_rows(df):
    return df[~df["tax_id"].isin(EMU_SPECIAL_ROWS)]


def available_ranks(emu):
    """Rangs présents dans toutes les tables (les bases n'ont pas toutes les mêmes)."""
    present = None
    for tables in emu.values():
        for df in tables.values():
            cols = {r for r in RANKS if r in df.columns}
            present = cols if present is None else present & cols
    return [r for r in RANKS if present and r in present]


def abundance_matrix(tables, rank):
    """Matrice taxons × échantillons (en %) à un rang donné."""
    series = {}
    for sample, df in tables.items():
        taxa = taxa_rows(df)
        if rank in taxa.columns:
            labels = taxa[rank].fillna("").astype(str).str.strip().replace("", UNKNOWN_LABEL)
        else:
            labels = taxa["tax_id"]
        series[sample] = taxa.groupby(labels)["abundance"].sum() * 100
    matrix = pd.DataFrame(series).fillna(0.0)
    return matrix[[s for s in tables]] if not matrix.empty else matrix


def classification_summary(tables):
    """Répartition des reads par catégorie Emu (en nombre de reads estimé)."""
    rows = {}
    for sample, df in tables.items():
        if "estimated counts" not in df.columns:
            continue
        special = df.set_index("tax_id")["estimated counts"]
        rows[sample] = {
            "Classés": taxa_rows(df)["estimated counts"].sum(),
            "Alignés non classés": special.get("mapped_unclassified", 0.0),
            "Alignés filtrés": special.get("mapped_filtered", 0.0),
            "Non alignés": special.get("unmapped", 0.0),
        }
    summary = pd.DataFrame.from_dict(rows, orient="index")
    if not summary.empty:
        summary["total"] = summary.sum(axis=1)
    return summary


def alpha_diversity(tables, rank):
    rows = {}
    matrix = abundance_matrix(tables, rank).drop(index=UNKNOWN_LABEL, errors="ignore")
    for sample, matrix_col in matrix.items():
        p = matrix_col[matrix_col > 0] / 100
        p = p / p.sum() if p.sum() > 0 else p
        rows[sample] = {
            "richness": int((p > 0).sum()),
            "shannon": float(-(p * np.log(p)).sum()) if len(p) else 0.0,
            "inv_simpson": float(1 / (p ** 2).sum()) if len(p) else 0.0,
        }
    return pd.DataFrame.from_dict(rows, orient="index")


def bray_curtis(a, b):
    """Dissimilarité de Bray-Curtis entre deux profils (en %)."""
    a, b = a.align(b, fill_value=0.0)
    total = a.sum() + b.sum()
    return float((a - b).abs().sum() / total) if total > 0 else np.nan


def top_taxa(matrices, top_n):
    """Taxons les plus abondants (moyenne sur tous les échantillons et bases)."""
    combined = pd.concat([m.mean(axis=1) for m in matrices if not m.empty], axis=1).fillna(0)
    ranked = combined.mean(axis=1).sort_values(ascending=False)
    ranked = ranked[ranked.index != UNKNOWN_LABEL]
    return list(ranked.index[:top_n])


def color_map(taxa):
    colors = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(taxa)}
    colors[OTHER_LABEL] = OTHER_COLOR
    colors[UNKNOWN_LABEL] = UNKNOWN_COLOR
    return colors


def collapse(matrix, keep):
    """Garde les taxons de 'keep' (+ Non assigné), regroupe le reste en 'Autres'."""
    keep = [t for t in keep if t in matrix.index]
    if UNKNOWN_LABEL in matrix.index:
        keep.append(UNKNOWN_LABEL)
    out = matrix.loc[keep].copy()
    rest = matrix.drop(index=keep).sum(axis=0)
    if (rest > 1e-9).any():
        out.loc[OTHER_LABEL] = rest
    return out


# ─────────────────────────────────────────────────────────────────────────────
# FIGURES
# ─────────────────────────────────────────────────────────────────────────────
def base_layout(fig, title, ytitle, height=480):
    fig.update_layout(
        title=title, height=height, template="plotly_white",
        yaxis_title=ytitle, legend=dict(traceorder="normal", font=dict(size=11)),
        margin=dict(l=60, r=20, t=60, b=100), hovermode="closest",
    )
    fig.update_xaxes(tickangle=-45)
    return fig


def fig_qc_reads(qc):
    fig = go.Figure()
    samples = list(qc.index)
    if "raw_number_of_reads" in qc:
        fig.add_bar(name="Reads bruts", x=samples, y=qc["raw_number_of_reads"], marker_color="#9EB9F3")
    if "kept_reads" in qc:
        fig.add_bar(name="Reads après chopper", x=samples, y=qc["kept_reads"], marker_color="#1F77B4")
    fig.update_layout(barmode="group")
    return base_layout(fig, "Nombre de reads avant / après filtrage", "Reads")


def fig_classification(summary, db):
    fig = go.Figure()
    pct = summary.drop(columns="total").div(summary["total"].replace(0, np.nan), axis=0) * 100
    for category in CLASSIF_COLORS:
        if category not in pct:
            continue
        fig.add_bar(
            name=category, x=list(pct.index), y=pct[category],
            marker_color=CLASSIF_COLORS[category],
            customdata=summary[category],
            hovertemplate="%{x}<br>" + category + " : %{y:.1f} % (%{customdata:.0f} reads)<extra></extra>",
        )
    fig.update_layout(barmode="stack", yaxis_range=[0, 100])
    return base_layout(fig, f"Devenir des reads — base « {db} »", "% des reads filtrés")


def fig_composition(matrix, colors, db, rank):
    fig = go.Figure()
    for taxon in matrix.index:
        fig.add_bar(
            name=taxon, x=list(matrix.columns), y=matrix.loc[taxon],
            marker_color=colors.get(taxon, OTHER_COLOR),
            hovertemplate="%{x}<br>" + esc(taxon) + " : %{y:.2f} %<extra></extra>",
        )
    fig.update_layout(barmode="stack", yaxis_range=[0, 100.5])
    return base_layout(fig, f"Composition — {RANK_LABELS.get(rank, rank)} — base « {db} »",
                       "Abondance relative (%)", height=560)


def fig_alpha(alpha_by_db, metric, label):
    fig = go.Figure()
    for i, (db, alpha) in enumerate(alpha_by_db.items()):
        if alpha.empty:
            continue
        fig.add_bar(name=db, x=list(alpha.index), y=alpha[metric], marker_color=PALETTE[i % len(PALETTE)])
    fig.update_layout(barmode="group")
    return base_layout(fig, label, label, height=420)


def fig_db_heatmap(matrices, taxa, rank):
    """Heatmap taxons × (échantillon | base) pour comparer les bases côte à côte."""
    samples = sorted({s for m in matrices.values() for s in m.columns}, key=natural_key)
    columns, labels = [], []
    for sample in samples:
        for db, matrix in matrices.items():
            col = matrix[sample] if sample in matrix.columns else pd.Series(dtype=float)
            columns.append(col.reindex(taxa).fillna(0.0))
            labels.append(f"{sample} | {db}")
    z = pd.concat(columns, axis=1).values if columns else np.zeros((0, 0))
    fig = go.Figure(go.Heatmap(
        z=z, x=labels, y=taxa, colorscale="YlGnBu", zmin=0,
        colorbar=dict(title="%"),
        hovertemplate="%{x}<br>%{y} : %{z:.2f} %<extra></extra>",
    ))
    fig.update_yaxes(autorange="reversed")
    height = max(420, 22 * len(taxa) + 200)
    return base_layout(fig, f"Comparaison des bases — {RANK_LABELS.get(rank, rank)} (top {len(taxa)})",
                       "", height=height)


def fig_bray_curtis(bc):
    fig = go.Figure()
    for i, pair in enumerate(bc.columns):
        fig.add_bar(name=pair, x=list(bc.index), y=bc[pair], marker_color=PALETTE[i % len(PALETTE)])
    fig.update_layout(barmode="group", yaxis_range=[0, 1])
    return base_layout(fig, "Dissimilarité de Bray-Curtis entre bases (0 = profils identiques)",
                       "Bray-Curtis", height=420)


# ─────────────────────────────────────────────────────────────────────────────
# HTML
# ─────────────────────────────────────────────────────────────────────────────
CSS = """
:root { --accent:#1F4E79; --bg:#F7F9FC; --border:#DDE3EA; --text:#1D2733; --muted:#5B6B7C; }
* { box-sizing: border-box; }
body { margin:0; font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       background: var(--bg); color: var(--text); }
header { background: var(--accent); color:#fff; padding: 22px 32px; }
header h1 { margin:0 0 4px 0; font-size: 1.6em; }
header p { margin:0; opacity:.85; }
nav { background:#fff; border-bottom:1px solid var(--border); padding: 8px 32px; position: sticky; top:0; z-index:10; }
nav a { margin-right: 18px; color: var(--accent); text-decoration:none; font-weight:600; font-size:.92em; }
main { max-width: 1400px; margin: 0 auto; padding: 16px 32px 48px; }
section { background:#fff; border:1px solid var(--border); border-radius:8px; padding: 18px 22px; margin: 18px 0; }
h2 { margin-top:0; color: var(--accent); }
.cards { display:flex; flex-wrap:wrap; gap:12px; margin: 8px 0 4px; }
.card { flex: 1 1 160px; background: var(--bg); border:1px solid var(--border); border-radius:8px; padding:12px 14px; }
.card .v { font-size:1.5em; font-weight:700; color: var(--accent); }
.card .l { font-size:.85em; color: var(--muted); }
table { border-collapse: collapse; width:100%; font-size:.88em; margin: 8px 0; }
th, td { border-bottom:1px solid var(--border); padding: 6px 8px; text-align:right; }
th { background: var(--bg); position: sticky; top: 0; }
th:first-child, td:first-child { text-align:left; }
td.txt { text-align:left; }
tr.excluded td { color:#A33; }
.tablewrap { max-height: 520px; overflow:auto; }
.note { color: var(--muted); font-size:.88em; }
.controls { margin: 6px 0 10px; }
.controls label { margin-right: 6px; font-weight:600; }
.controls select { margin-right: 18px; padding: 3px 6px; }
.panel { display:none; }
.panel.active { display:block; }
.badge { display:inline-block; padding: 1px 7px; border-radius: 10px; font-size:.8em; }
.badge.ok { background:#E3F4E3; color:#1E6B1E; }
.badge.ko { background:#FBE3E3; color:#A33; }
code { background: var(--bg); padding: 1px 4px; border-radius:3px; }
@media print { nav { display:none; } .panel { display:block !important; } section { break-inside: avoid; } }
"""

JS = """
function showPanel(group) {
  var selects = document.querySelectorAll('select[data-group="' + group + '"]');
  var key = Array.prototype.map.call(selects, function (s) { return s.value; }).join('__');
  document.querySelectorAll('.panel[data-group="' + group + '"]').forEach(function (p) {
    var active = p.getAttribute('data-key') === key;
    p.classList.toggle('active', active);
    if (active) {
      p.querySelectorAll('.js-plotly-plot').forEach(function (g) { Plotly.Plots.resize(g); });
    }
  });
}
document.addEventListener('DOMContentLoaded', function () {
  var groups = {};
  document.querySelectorAll('select[data-group]').forEach(function (s) {
    groups[s.getAttribute('data-group')] = true;
    s.addEventListener('change', function () { showPanel(s.getAttribute('data-group')); });
  });
  Object.keys(groups).forEach(showPanel);
});
"""


def select_html(group, label, options, selected=None):
    opts = "".join(
        f'<option value="{esc(slug(v))}"{" selected" if v == selected else ""}>{esc(text)}</option>'
        for v, text in options
    )
    return f'<label>{esc(label)}</label><select data-group="{group}">{opts}</select>'


def panel_html(group, key_parts, content):
    key = "__".join(slug(k) for k in key_parts)
    return f'<div class="panel" data-group="{group}" data-key="{esc(key)}">{content}</div>'


def table_html(headers, rows, row_classes=None):
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for i, row in enumerate(rows):
        cls = f' class="{row_classes[i]}"' if row_classes and row_classes[i] else ""
        cells = "".join(f'<td class="txt">{c}</td>' if j == 0 else f"<td>{c}</td>" for j, c in enumerate(row))
        body.append(f"<tr{cls}>{cells}</tr>")
    return f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def section_summary(qc, emu, run_info):
    n_samples = len(qc.index) if not qc.empty else 0
    n_raw = qc["raw_number_of_reads"].sum() if "raw_number_of_reads" in qc else np.nan
    n_kept = qc["kept_reads"].sum() if "kept_reads" in qc else np.nan
    n_classified = len({s for tables in emu.values() for s in tables})
    cards = [
        (fmt_int(n_samples), "échantillons"),
        (fmt_int(n_classified), "échantillons classés par Emu"),
        (fmt_int(n_raw), "reads bruts"),
        (fmt_int(n_kept), "reads après filtrage"),
        (str(len(emu)), "base(s) Emu : " + ", ".join(esc(d) for d in emu) if emu else "aucune base"),
    ]
    cards_html = "".join(f'<div class="card"><div class="v">{v}</div><div class="l">{l}</div></div>'
                         for v, l in cards)
    labels = {
        "pipeline": "Pipeline", "date": "Date", "fastq_dir": "Dossier FASTQ",
        "sample_sheet": "Sample sheet", "min_quality": "Qualité min. (chopper)",
        "min_length": "Longueur min. (pb)", "max_length": "Longueur max. (pb)",
        "chopper_args": "Options chopper", "min_reads": "Reads min. pour Emu",
        "emu_dbs": "Bases Emu", "emu_min_abundance": "Abondance min. Emu",
        "emu_N": "Alignements secondaires (--N)", "emu_args": "Options Emu",
        "resultsdir": "Dossier des résultats",
    }
    params_rows = [[esc(labels.get(k, k)), f"<code>{esc(v)}</code>"] for k, v in run_info.items()]
    return f"""
<section id="resume">
  <h2>Résumé</h2>
  <div class="cards">{cards_html}</div>
  <details><summary>Paramètres du run</summary>{table_html(["Paramètre", "Valeur"], params_rows)}</details>
</section>"""


def section_qc(qc, emu, min_reads):
    if qc.empty:
        return '<section id="qc"><h2>Contrôle qualité</h2><p class="note">Aucune statistique NanoPlot disponible.</p></section>'
    classified = {s for tables in emu.values() for s in tables}
    rows, classes = [], []
    for sample, r in qc.iterrows():
        kept = r.get("kept_reads", np.nan)
        excluded = not (isinstance(kept, (int, float)) and not np.isnan(kept) and kept >= min_reads)
        status = ('<span class="badge ko">exclu (&lt; {} reads)</span>'.format(min_reads) if excluded
                  else '<span class="badge ok">classé</span>' if sample in classified
                  else '<span class="badge ko">non classé</span>')
        rows.append([
            esc(sample),
            fmt_int(r.get("raw_number_of_reads")), fmt_int(kept), fmt_float(r.get("pct_kept")),
            fmt_int(r.get("raw_median_read_length")), fmt_int(r.get("filtered_median_read_length")),
            fmt_float(r.get("raw_median_qual")), fmt_float(r.get("filtered_median_qual")),
            fmt_int(r.get("filtered_n50")), status,
        ])
        classes.append("excluded" if excluded else "")
    headers = ["Échantillon", "Reads bruts", "Reads filtrés", "% conservés",
               "Long. médiane brute", "Long. médiane filtrée", "Q médiane brute",
               "Q médiane filtrée", "N50 filtré", "Statut"]
    return f"""
<section id="qc">
  <h2>Contrôle qualité des reads (NanoPlot + chopper)</h2>
  {fig_html(fig_qc_reads(qc), "fig_qc_reads")}
  {table_html(headers, rows, classes)}
  <p class="note">Rapports NanoPlot détaillés : <code>results/nanoplot/&lt;échantillon&gt;_raw</code> et
  <code>_filtered</code> (fichier <code>NanoPlot-report.html</code>). Les échantillons avec moins de
  {min_reads} reads après filtrage ne sont pas classés par Emu.</p>
</section>"""


def section_classification(emu):
    if not emu:
        return ""
    dbs = list(emu)
    panels, table_rows = [], []
    for db in dbs:
        summary = classification_summary(emu[db])
        if summary.empty:
            continue
        panels.append(panel_html("classif", [db], fig_html(fig_classification(summary, db), f"fig_classif_{slug(db)}")))
        for sample, r in summary.iterrows():
            pct = 100 * r["Classés"] / r["total"] if r["total"] else np.nan
            table_rows.append([esc(sample), esc(db), fmt_int(r["total"]), fmt_int(r["Classés"]),
                               fmt_float(pct), fmt_int(r["Alignés non classés"]),
                               fmt_int(r["Alignés filtrés"]), fmt_int(r["Non alignés"])])
    table_rows.sort(key=lambda row: natural_key(row[0]))
    headers = ["Échantillon", "Base", "Reads", "Classés", "% classés",
               "Alignés non classés", "Alignés filtrés", "Non alignés"]
    return f"""
<section id="classification">
  <h2>Classification Emu</h2>
  <div class="controls">{select_html("classif", "Base :", [(d, d) for d in dbs])}</div>
  {"".join(panels)}
  {table_html(headers, table_rows)}
  <p class="note"><b>Non alignés</b> : aucun alignement minimap2 sur la base.
  <b>Alignés filtrés</b> : alignements rejetés par les seuils de longueur / identité.
  <b>Alignés non classés</b> : reads alignés uniquement sur des espèces jugées absentes par Emu
  (trop faible abondance). Comptages estimés par Emu.</p>
</section>"""


def section_composition(emu, ranks, top_n):
    if not emu or not ranks:
        return ""
    panels = []
    for rank in ranks:
        matrices = {db: abundance_matrix(tables, rank) for db, tables in emu.items()}
        colors = color_map(top_taxa(list(matrices.values()), top_n * 2))
        for db, matrix in matrices.items():
            keep = top_taxa([matrix], top_n)
            collapsed = collapse(matrix, keep)
            panels.append(panel_html("compo", [db, rank],
                                     fig_html(fig_composition(collapsed, colors, db, rank),
                                              f"fig_compo_{slug(db)}_{rank}")))
    default_rank = "genus" if "genus" in ranks else ranks[0]
    controls = (select_html("compo", "Base :", [(d, d) for d in emu]) +
                select_html("compo", "Rang :", [(r, RANK_LABELS[r]) for r in ranks], default_rank))
    return f"""
<section id="composition">
  <h2>Composition taxonomique</h2>
  <div class="controls">{controls}</div>
  {"".join(panels)}
  <p class="note">Les {top_n} taxons les plus abondants sont affichés, les autres sont regroupés dans
  « {OTHER_LABEL} ». Abondances relatives calculées sur les reads classés. Un même taxon garde la même
  couleur d'une base à l'autre. Tables complètes : <code>results/emu/&lt;base&gt;/</code>.</p>
</section>"""


def section_alpha(emu, ranks):
    if not emu or not ranks:
        return ""
    panels, rows = [], []
    for rank in ranks:
        alpha_by_db = {db: alpha_diversity(tables, rank) for db, tables in emu.items()}
        figs = "".join(fig_html(fig_alpha(alpha_by_db, metric, label), f"fig_alpha_{rank}_{metric}")
                       for metric, label in [("shannon", "Indice de Shannon"),
                                             ("richness", "Richesse (nb de taxons)")])
        panels.append(panel_html("alpha", [rank], figs))
        if rank == ranks[0]:
            for db, alpha in alpha_by_db.items():
                for sample, r in alpha.iterrows():
                    rows.append([esc(sample), esc(db), fmt_int(r["richness"]),
                                 fmt_float(r["shannon"], 2), fmt_float(r["inv_simpson"], 2)])
    rows.sort(key=lambda row: natural_key(row[0]))
    return f"""
<section id="diversite">
  <h2>Diversité alpha</h2>
  <div class="controls">{select_html("alpha", "Rang :", [(r, RANK_LABELS[r]) for r in ranks])}</div>
  {"".join(panels)}
  <h3>Au rang {RANK_LABELS[ranks[0]].lower()}</h3>
  {table_html(["Échantillon", "Base", "Richesse", "Shannon", "Simpson inverse"], rows)}
  <p class="note">La richesse dépend du seuil <code>--emu_min_abundance</code> et de la profondeur de
  séquençage : la comparer entre échantillons de profondeurs proches.</p>
</section>"""


def section_db_comparison(emu, ranks, top_n):
    if len(emu) < 2 or not ranks:
        return """
<section id="comparaison">
  <h2>Comparaison des bases</h2>
  <p class="note">Une seule base utilisée. Pour comparer, relancer avec plusieurs bases,
  par exemple <code>--emu_dbs emu,silva,rdp</code>.</p>
</section>"""
    dbs = list(emu)
    pairs = [(a, b) for i, a in enumerate(dbs) for b in dbs[i + 1:]]
    panels = []
    for rank in ranks:
        matrices = {db: abundance_matrix(tables, rank) for db, tables in emu.items()}
        taxa = top_taxa(list(matrices.values()), top_n * 2)
        samples = sorted({s for m in matrices.values() for s in m.columns}, key=natural_key)

        # Bray-Curtis par échantillon et par paire de bases
        bc = pd.DataFrame(index=samples, columns=[f"{a} vs {b}" for a, b in pairs], dtype=float)
        for a, b in pairs:
            for s in samples:
                if s in matrices[a].columns and s in matrices[b].columns:
                    bc.loc[s, f"{a} vs {b}"] = bray_curtis(matrices[a][s], matrices[b][s])

        # Taxon dominant par échantillon et par base
        dom_rows = []
        for s in samples:
            row = [esc(s)]
            for db in dbs:
                m = matrices[db]
                if s in m.columns and m[s].sum() > 0:
                    col = m[s].drop(index=UNKNOWN_LABEL, errors="ignore")
                    top = col.idxmax() if not col.empty else UNKNOWN_LABEL
                    row.append(f"{esc(top)} ({fmt_float(m.loc[top, s])} %)")
                else:
                    row.append("-")
            dom_rows.append(row)

        bc_rows = [[esc(s)] + [fmt_float(v, 3) for v in bc.loc[s]] for s in samples]
        content = (
            fig_html(fig_db_heatmap(matrices, taxa, rank), f"fig_dbheat_{rank}") +
            fig_html(fig_bray_curtis(bc), f"fig_dbbc_{rank}") +
            "<h3>Bray-Curtis entre bases</h3>" + table_html(["Échantillon"] + list(bc.columns), bc_rows) +
            "<h3>Taxon dominant</h3>" + table_html(["Échantillon"] + dbs, dom_rows)
        )
        panels.append(panel_html("dbcomp", [rank], content))
    default_rank = "genus" if "genus" in ranks else ranks[0]
    return f"""
<section id="comparaison">
  <h2>Comparaison des bases</h2>
  <div class="controls">{select_html("dbcomp", "Rang :", [(r, RANK_LABELS[r]) for r in ranks], default_rank)}</div>
  {"".join(panels)}
  <p class="note">Les noms de taxons sont comparés tels quels : une différence de nomenclature entre bases
  (ex. <i>Escherichia-Shigella</i> dans SILVA, <i>Escherichia</i> dans NCBI) apparaît comme une
  différence de composition. La comparaison au rang genre est en général la plus pertinente.</p>
</section>"""


def build_report(qc, emu, run_info, title, top_n, min_reads):
    ranks = available_ranks(emu)
    sections = [
        section_summary(qc, emu, run_info),
        section_qc(qc, emu, min_reads),
        section_classification(emu),
        section_composition(emu, ranks, top_n),
        section_alpha(emu, ranks),
        section_db_comparison(emu, ranks, top_n),
    ]
    nav = [("resume", "Résumé"), ("qc", "QC"), ("classification", "Classification"),
           ("composition", "Composition"), ("diversite", "Diversité"), ("comparaison", "Comparaison des bases")]
    nav_html = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in nav)
    subtitle = f"{esc(run_info.get('pipeline', '16StJo'))} — {esc(run_info.get('date', ''))}"
    return f"""<!DOCTYPE html>
<html lang="fr">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<style>{CSS}</style>
<script type="text/javascript">{get_plotlyjs()}</script>
</head>
<body>
<header><h1>{esc(title)}</h1><p>{subtitle}</p></header>
<nav>{nav_html}</nav>
<main>
{"".join(sections)}
</main>
<script>{JS}</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Rapport HTML du pipeline 16StJo (Emu)")
    parser.add_argument("--nanostats_dir", required=True)
    parser.add_argument("--read_counts", default=None)
    parser.add_argument("--emu_dir", required=True)
    parser.add_argument("--run_info", default=None)
    parser.add_argument("--title", default="Analyse 16S Nanopore - Emu")
    parser.add_argument("--top_n", type=int, default=15)
    parser.add_argument("--min_reads", type=int, default=50)
    parser.add_argument("--output", default="16StJo_report.html")
    args = parser.parse_args()

    run_info = {}
    if args.run_info and os.path.exists(args.run_info):
        with open(args.run_info) as handle:
            run_info = json.load(handle)

    qc = load_qc(args.nanostats_dir, args.read_counts)
    emu = load_emu(args.emu_dir)
    report = build_report(qc, emu, run_info, args.title, args.top_n, args.min_reads)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(report)
    print(f"Rapport écrit : {args.output} ({len(qc)} échantillons, {len(emu)} base(s))")


if __name__ == "__main__":
    main()
