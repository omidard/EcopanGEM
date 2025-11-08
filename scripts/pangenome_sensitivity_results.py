#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pangenome_sensitivity_car_dotplots.py

What you get (publication-ready; PNG + SVG):
(A) clusters_vs_cdhit_threshold.(png|svg) — Bar plot of total cluster counts vs CD-HIT identity threshold.
(B) CAR_sensitivity_cdhit80.(png|svg) — Three *dot plots* (Core / Accessory / Rare) for CD-HIT 80%:
    • Left column: vary CORE threshold 90→100 (RARE fixed at 6.8%)
    • Right column: vary RARE threshold 10→1 (CORE fixed at 96.7%)
    Each panel shows dots (with a faint guide line), carefully formatted ticks/labels/legends.

Plus three standalone category figures (dot plots with the two sweeps per category):
    CAR_core_sensitivity_cdhit80.(png|svg)
    CAR_accessory_sensitivity_cdhit80.(png|svg)
    CAR_rare_sensitivity_cdhit80.(png|svg)

Key fixes:
- **Strict genome-ID parsing** (no more NP_418157.1 “fake genomes”).
- **TOTAL_80** taken from the presence/absence matrix columns when available
  (tries `presence_absence_matrix_80.csv`, then `presence_absence_matrix.csv`).
- Multiprocessing across **all cores** + `tqdm` progress bars.
- NA-safe plotting (drops thresholds with missing JSONs for panel A).

Run:
    python pangenome_sensitivity_car_dotplots.py
Edit BASE if needed.
"""

from __future__ import annotations
import json
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm

# -------------------------
# CONFIG
# -------------------------
BASE = Path("/home/omidard/ecoli_revision/pangenome_s")   # where cluster_to_locus_*.json + presence_absence_matrix_*.csv live
THRESHOLDS = [65, 70, 75, 80, 85, 90, 95]

# Benchmark CAR thresholds (from your inflection method)
CORE_BENCH = 96.7
RARE_BENCH = 6.8

# Sensitivity sweeps
CORE_SWEEP = list(range(90, 101))      # 90 → 100 inclusive
RARE_SWEEP = list(range(10, 0, -1))    # 10 → 1 inclusive (descending on purpose)

OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- outputs ----
PANEL_A_PNG = OUT_DIR / "clusters_vs_cdhit_threshold.png"
PANEL_A_SVG = OUT_DIR / "clusters_vs_cdhit_threshold.svg"

PANEL_B_PNG = OUT_DIR / "CAR_sensitivity_cdhit80.png"
PANEL_B_SVG = OUT_DIR / "CAR_sensitivity_cdhit80.svg"

PANEL_B_CORE_PNG = OUT_DIR / "CAR_core_sensitivity_cdhit80.png"
PANEL_B_CORE_SVG = OUT_DIR / "CAR_core_sensitivity_cdhit80.svg"
PANEL_B_ACC_PNG  = OUT_DIR / "CAR_accessory_sensitivity_cdhit80.png"
PANEL_B_ACC_SVG  = OUT_DIR / "CAR_accessory_sensitivity_cdhit80.svg"
PANEL_B_RARE_PNG = OUT_DIR / "CAR_rare_sensitivity_cdhit80.png"
PANEL_B_RARE_SVG = OUT_DIR / "CAR_rare_sensitivity_cdhit80.svg"

# -------------------------
# Helpers (STRICT parsing)
# -------------------------
GCF_GCA_RE = re.compile(r"(GCA|GCF)_[0-9]+(?:\.[0-9]+)?")

def json_candidates(thr: int) -> List[Path]:
    cands = [BASE / f"cluster_to_locus_{thr}.json"]
    if thr == 80:
        cands.append(BASE / "cluster_to_locus.json")  # legacy fallback
    return cands

def find_json_for_threshold(thr: int) -> Path:
    for cand in json_candidates(thr):
        if cand.is_file():
            return cand
    raise FileNotFoundError(f"No JSON for threshold {thr}. Tried: {[str(c) for c in json_candidates(thr)]}")

def load_clusters_as_genome_sets(json_path: Path) -> Dict[str, Set[str]]:
    """
    Normalize JSON to: { cluster_id : set(genome_ids) }

    Allowed sources ONLY:
      - dict members with 'genome_id'
      - string members like '...|GENOME|...'
      - last-resort: a GCF_/GCA_ token anywhere in the record

    We NEVER treat generic 'digits.digits' (e.g., NP_418157.1) as genomes.
    """
    with json_path.open("r") as f:
        raw = json.load(f)

    clusters: Dict[str, Set[str]] = {}
    for cid, members in raw.items():
        gset: Set[str] = set()
        for m in (members or []):
            if isinstance(m, dict):
                gid = str(m.get("genome_id", "")).strip()
                if gid:
                    gset.add(gid)
                else:
                    token = None
                    for v in m.values():
                        if isinstance(v, str):
                            mo = GCF_GCA_RE.search(v)
                            if mo:
                                token = mo.group(0)
                                break
                    if token:
                        gset.add(token)
            elif isinstance(m, str):
                parts = m.strip().lstrip(">").split("|")
                if len(parts) >= 2 and parts[1].strip():
                    gset.add(parts[1].strip())
                else:
                    mo = GCF_GCA_RE.search(m)
                    if mo:
                        gset.add(mo.group(0))
                # else: skip silently (do NOT guess)
        if gset:
            clusters[str(cid)] = gset
    return clusters

def classify_counts(
    clusters: Dict[str, Set[str]],
    core_thr: float,
    rare_thr: float,
    total_genomes: int
) -> Tuple[int, int, int]:
    """
    Return (core_count, accessory_count, rare_count), thresholds in percentage.
      core if frac > core_thr
      rare if frac < rare_thr
      accessory otherwise
    """
    core_n = acc_n = rare_n = 0
    for gset in clusters.values():
        frac = (len(gset) / total_genomes) * 100.0 if total_genomes else 0.0
        if frac > core_thr:
            core_n += 1
        elif frac < rare_thr:
            rare_n += 1
        else:
            acc_n += 1
    return core_n, acc_n, rare_n

# -------------------------
# MP Tasks
# -------------------------
def task_panel_a(thr: int) -> Tuple[int, int]:
    try:
        jp = find_json_for_threshold(thr)
        clus = load_clusters_as_genome_sets(jp)
        return (thr, len(clus))
    except FileNotFoundError:
        return (thr, -1)

def task_core_sweep(cthr: int, clusters_80: Dict[str, Set[str]], total_80: int) -> Tuple[int, int, int, int]:
    core_n, acc_n, rare_n = classify_counts(clusters_80, core_thr=cthr, rare_thr=RARE_BENCH, total_genomes=total_80)
    return (cthr, core_n, acc_n, rare_n)

def task_rare_sweep(rthr: int, clusters_80: Dict[str, Set[str]], total_80: int) -> Tuple[int, int, int, int]:
    core_n, acc_n, rare_n = classify_counts(clusters_80, core_thr=CORE_BENCH, rare_thr=rthr, total_genomes=total_80)
    return (rthr, core_n, acc_n, rare_n)

# -------------------------
# Plot helpers (million-dollar polish)
# -------------------------
def set_pubstyle():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "legend.fontsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

# color palette (color-blind friendly)
COL_RARE = "#1f78b4"      # blue
COL_ACC  = "#ff7f00"      # orange
COL_CORE = "#6a3d9a"      # purple
GUIDE_COLOR = "#666666"

def dotplot(ax, x, y, color, label, marker="o"):
    ax.plot(x, y, marker=marker, linestyle="None", ms=6, mec="black", mfc=color, label=label)
    ax.plot(x, y, linestyle="-", linewidth=1, alpha=0.4, color=GUIDE_COLOR)  # faint guide line

def fmt_axis_common(ax, xlabel, ylabel="Cluster count"):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_formatter(mtick.StrMethodFormatter('{x:,.0f}'))
    ax.grid(axis="y", linestyle=":", linewidth=0.7, alpha=0.6)

# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    nproc = cpu_count()
    print(f"Using {nproc} CPU cores.")

    # ------------------ Panel (A): clusters vs threshold ------------------
    with Pool(processes=nproc) as pool:
        results_a = list(tqdm(pool.imap_unordered(task_panel_a, THRESHOLDS),
                              total=len(THRESHOLDS), desc="Panel (A) counting clusters"))
    df_a = pd.DataFrame(results_a, columns=["threshold", "clusters"]).sort_values("threshold")
    df_a["clusters"] = df_a["clusters"].replace({-1: pd.NA})
    df_a_plot = df_a.dropna(subset=["clusters"]).copy()
    df_a_plot["threshold"] = pd.to_numeric(df_a_plot["threshold"], errors="coerce").astype(int)
    df_a_plot["clusters"]  = pd.to_numeric(df_a_plot["clusters"],  errors="coerce").astype(int)

    set_pubstyle()
    figA, axA = plt.subplots(figsize=(7.6, 4.2))
    xA = df_a_plot["threshold"].tolist()
    yA = df_a_plot["clusters"].tolist()
    bars = axA.bar(xA, yA, width=4, edgecolor="black", facecolor="#9ecae1")
    axA.set_xlabel("CD-HIT identity threshold (%)")
    axA.set_ylabel("Number of protein clusters")
    axA.set_title("Total clusters vs CD-HIT identity")
    axA.set_xticks(xA)
    axA.set_xticklabels([str(t) for t in xA])
    axA.yaxis.set_major_formatter(mtick.StrMethodFormatter('{x:,.0f}'))
    for b in bars:
        h = b.get_height()
        axA.annotate(f"{int(h):,}", xy=(b.get_x()+b.get_width()/2, h),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", fontsize=9)
    axA.text(-0.12, 1.06, "(A)", transform=axA.transAxes, fontweight="bold", fontsize=16)
    figA.tight_layout()
    figA.savefig(PANEL_A_PNG, bbox_inches="tight")
    figA.savefig(PANEL_A_SVG, bbox_inches="tight")
    plt.close(figA)

    # ------------------ Panel (B): CAR @ 80% as DOT PLOTS -----------------
    # Load clusters@80
    json_80 = find_json_for_threshold(80)
    clusters_80 = load_clusters_as_genome_sets(json_80)

    # Decide TOTAL_80 from PA matrix columns if available
    pa_80 = BASE / "presence_absence_matrix_80.csv"
    pa_any = BASE / "presence_absence_matrix.csv"
    if pa_80.is_file():
        TOTAL_80 = len(pd.read_csv(pa_80, nrows=0).columns)
    elif pa_any.is_file():
        TOTAL_80 = len(pd.read_csv(pa_any, nrows=0).columns)
    else:
        # fallback to JSON (unique genome IDs)
        TOTAL_80 = len(set().union(*clusters_80.values())) if clusters_80 else 0

    print(f"[INFO] Total genomes for CD-HIT 80% (denominator): {TOTAL_80}")

    # Sanity check: this should match your PA columns (should be 2377)
    if TOTAL_80 and TOTAL_80 != 2377:
        print(f"[WARN] TOTAL_80={TOTAL_80} != 2377. If PA exists, please point the path correctly.")

    core_b, acc_b, rare_b = classify_counts(clusters_80, CORE_BENCH, RARE_BENCH, TOTAL_80)

    # Compute sweeps in parallel
    with Pool(processes=nproc) as pool:
        core_iter = pool.imap_unordered(partial(task_core_sweep, clusters_80=clusters_80, total_80=TOTAL_80), CORE_SWEEP)
        results_core = list(tqdm(core_iter, total=len(CORE_SWEEP), desc="Panel (B) core sweep"))
    df_core = pd.DataFrame(results_core, columns=["thr", "core", "accessory", "rare"]).sort_values("thr")

    with Pool(processes=nproc) as pool:
        rare_iter = pool.imap_unordered(partial(task_rare_sweep, clusters_80=clusters_80, total_80=TOTAL_80), RARE_SWEEP)
        results_rare = list(tqdm(rare_iter, total=len(RARE_SWEEP), desc="Panel (B) rare sweep"))
    df_rare = pd.DataFrame(results_rare, columns=["thr", "core", "accessory", "rare"]).sort_values("thr", ascending=False)

    bench_txt = (
        f"Benchmark: core>{CORE_BENCH:.1f}%, rare<{RARE_BENCH:.1f}%  →  "
        f"core={core_b:,}, accessory={acc_b:,}, rare={rare_b:,}   "
        f"(Total genomes used: {TOTAL_80:,})"
    )

    # ---- B (compact multi-panel): 3 categories × 2 sweeps (dot plots) ----
    set_pubstyle()
    figB, axes = plt.subplots(3, 2, figsize=(11.5, 9.6), sharey="row")
    categories = [("Core clusters", "core", COL_CORE),
                  ("Accessory clusters", "accessory", COL_ACC),
                  ("Rare clusters", "rare", COL_RARE)]
    for row, (title, colname, color) in enumerate(categories):
        # Left: vary core (rare fixed)
        axL = axes[row, 0]
        dotplot(axL, df_core["thr"].tolist(), df_core[colname].tolist(), color=color, label=f"{title} (vary core)")
        fmt_axis_common(axL, "Core threshold (%)  [rare fixed at 6.8%]")
        axL.set_title(title)
        axL.set_xticks(CORE_SWEEP[::2])
        axL.axvline(CORE_BENCH, color=GUIDE_COLOR, linestyle=":", linewidth=1)
        if row == 0:
            axL.legend(loc="upper left", frameon=True)

        # Right: vary rare (core fixed), x shown as 10→1
        axR = axes[row, 1]
        dotplot(axR, df_rare["thr"].tolist(), df_rare[colname].tolist(), color=color, label=f"{title} (vary rare)", marker="^")
        fmt_axis_common(axR, "Rare threshold (%)  [core fixed at 96.7%]")
        axR.set_xticks(RARE_SWEEP[::2])
        if row == 0:
            axR.legend(loc="upper left", frameon=True)

    axes[0, 0].text(-0.18, 1.18, "(B)", transform=axes[0, 0].transAxes, fontsize=16, fontweight="bold")
    axes[-1, 0].text(0.01, -0.28, bench_txt, transform=axes[-1, 0].transAxes, fontsize=10, ha="left", va="top")
    figB.tight_layout()
    figB.savefig(PANEL_B_PNG, bbox_inches="tight")
    figB.savefig(PANEL_B_SVG, bbox_inches="tight")
    plt.close(figB)

    # ---- Standalone category figures (requested) ----
    def save_category_figure(name: str, color: str, y_core: List[int], y_rare: List[int], out_png: Path, out_svg: Path):
        set_pubstyle()
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
        # Left
        dotplot(axL, CORE_SWEEP, y_core, color=color, label=f"{name} (vary core)")
        fmt_axis_common(axL, "Core threshold (%)  [rare fixed at 6.8%]")
        axL.set_title(f"{name} vs core threshold (CD-HIT 80%)")
        axL.set_xticks(CORE_SWEEP[::2])
        axL.axvline(CORE_BENCH, color=GUIDE_COLOR, linestyle=":", linewidth=1)
        axL.legend(loc="upper left", frameon=True)
        axL.text(-0.17, 1.12, "(B)", transform=axL.transAxes, fontsize=16, fontweight="bold")
        # Right
        dotplot(axR, RARE_SWEEP, y_rare, color=color, label=f"{name} (vary rare)", marker="^")
        fmt_axis_common(axR, "Rare threshold (%)  [core fixed at 96.7%]")
        axR.set_title(f"{name} vs rare threshold (CD-HIT 80%)")
        axR.set_xticks(RARE_SWEEP[::2])
        axR.legend(loc="upper left", frameon=True)
        axL.text(0.01, -0.20, bench_txt, transform=axL.transAxes, fontsize=10, ha="left", va="top")
        fig.tight_layout()
        fig.savefig(out_png, bbox_inches="tight")
        fig.savefig(out_svg, bbox_inches="tight")
        plt.close(fig)

    save_category_figure(
        "Core clusters", COL_CORE,
        df_core["core"].tolist(), df_rare["core"].tolist(),
        PANEL_B_CORE_PNG, PANEL_B_CORE_SVG
    )
    save_category_figure(
        "Accessory clusters", COL_ACC,
        df_core["accessory"].tolist(), df_rare["accessory"].tolist(),
        PANEL_B_ACC_PNG, PANEL_B_ACC_SVG
    )
    save_category_figure(
        "Rare clusters", COL_RARE,
        df_core["rare"].tolist(), df_rare["rare"].tolist(),
        PANEL_B_RARE_PNG, PANEL_B_RARE_SVG
    )

    # ------------------ Save numeric tables for reproducibility ------------
    df_a.to_csv(OUT_DIR / "clusters_vs_cdhit_threshold.csv", index=False)
    df_core.to_csv(OUT_DIR / "CAR_sensitivity_core_sweep.csv", index=False)
    df_rare.to_csv(OUT_DIR / "CAR_sensitivity_rare_sweep.csv", index=False)

    print(f"[OK] Panel (A): {PANEL_A_PNG} | {PANEL_A_SVG}")
    print(f"[OK] Panel (B): {PANEL_B_PNG} | {PANEL_B_SVG}")
    print(f"[OK] Standalone: {PANEL_B_CORE_PNG}, {PANEL_B_ACC_PNG}, {PANEL_B_RARE_PNG}")
    print(f"[OK] Tables in: {OUT_DIR}")
