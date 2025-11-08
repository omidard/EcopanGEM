#!/usr/bin/env python3
"""
Filter genomes by quality using robust stats + hard biological guards.

INPUT
  - A QC table with columns:
      genome, num_contigs, N50, gc_percent, completeness, contamination,
      fasta_path, gbff_path
  - A BASE directory where each genome lives in its own subfolder named like the
    'genome' column (e.g., .../annotations/GCA_000000000.1_ASM...).

DECISION RULES (by default; all tunable via CLI):
  Hard guards (biologically motivated; applied first):
    - completeness >= 90
    - contamination <= 5

  Robust statistical guards (median/MAD-based; applied to what remains):
    - num_contigs: reject if robust z > 2.5  (too many contigs)
    - N50:        reject if robust z < -2.5  (too small N50)
    - gc_percent: reject if |robust z| > 3.5 (GC outlier)
    - completeness: reject if robust z < -2.5 (unusually low)
    - contamination: reject if robust z > 2.5 (unusually high)

  Composite rule:
    - A genome fails if it fails ANY of the above criteria.

ACTIONS:
  - Default is DRY RUN: prints what would be moved/deleted, changes nothing.
  - Use --apply to actually move/remove folders.
  - By default, failing genomes are moved to --quarantine (created if needed).
  - Use --delete to permanently delete (dangerous!). Mutually exclusive with --quarantine.

USAGE
  python filter_genomes.py \
    --qc /path/to/genome_qc_results.csv \
    --base /home/omidard/panGEM_pipeline/ESKAPE/Enterococcus_faecium/annotations \
    --apply \
    --quarantine /home/omidard/panGEM_pipeline/ESKAPE/Enterococcus_faecium/low_quality

You can tweak thresholds; run with -h to see options.
"""

import argparse
import os
import sys
import shutil
from typing import Tuple, Optional

import numpy as np
import pandas as pd


def robust_z(x: pd.Series) -> pd.Series:
    """Median/MAD robust z-score. Returns NaN where MAD==0 or x is NaN."""
    x = pd.to_numeric(x, errors="coerce")
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad == 0:
        return pd.Series(np.full(len(x), np.nan), index=x.index)
    # 1.4826 scales MAD to be comparable to std for normal data
    return (x - med) / (1.4826 * mad)


def must_have_columns(df: pd.DataFrame, cols) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns in QC table: {missing}")


def safe_genome_folder(base_dir: str, genome: str, gbff_path: Optional[str]) -> Optional[str]:
    """
    Resolve the on-disk folder for a genome.
    Prefer the parent of gbff_path if it lives under base_dir; otherwise assume base_dir/genome.
    Returns None if nothing sensible is found.
    """
    if gbff_path and isinstance(gbff_path, str) and os.path.isabs(gbff_path):
        parent = os.path.dirname(gbff_path)
        try:
            common = os.path.commonpath([os.path.realpath(parent), os.path.realpath(base_dir)])
            if common == os.path.realpath(base_dir):
                return parent
        except Exception:
            pass
    # Fallback to base_dir/genome
    candidate = os.path.join(base_dir, genome)
    if os.path.isdir(candidate):
        return candidate
    return None


def explain_fail(row: pd.Series, cfg) -> Tuple[bool, list]:
    """
    Return (failed, reasons[]) for a genome row.
    Uses both hard guards and robust z-score guards.
    """
    reasons = []

    # Hard guards
    comp = pd.to_numeric(row.get("completeness"), errors="coerce")
    cont = pd.to_numeric(row.get("contamination"), errors="coerce")

    if np.isfinite(cfg.comp_min) and (not np.isfinite(comp) or comp < cfg.comp_min):
        reasons.append(f"completeness<{cfg.comp_min} ({comp})")

    if np.isfinite(cfg.cont_max) and (not np.isfinite(cont) or cont > cfg.cont_max):
        reasons.append(f"contamination>{cfg.cont_max} ({cont})")

    # Robust z (these already computed and present in row)
    if pd.notna(row.get("num_contigs_rz")) and row["num_contigs_rz"] > cfg.num_contigs_rz_max:
        reasons.append(f"num_contigs_rz>{cfg.num_contigs_rz_max:.2f} ({row['num_contigs_rz']:.2f})")
    if pd.notna(row.get("N50_rz")) and row["N50_rz"] < -cfg.N50_rz_max_abs:
        reasons.append(f"N50_rz<-{cfg.N50_rz_max_abs:.2f} ({row['N50_rz']:.2f})")
    if pd.notna(row.get("gc_percent_rz")) and abs(row["gc_percent_rz"]) > cfg.gc_rz_max_abs:
        reasons.append(f"|gc_percent_rz|>{cfg.gc_rz_max_abs:.2f} ({row['gc_percent_rz']:.2f})")
    if pd.notna(row.get("completeness_rz")) and row["completeness_rz"] < -cfg.comp_rz_max_abs:
        reasons.append(f"completeness_rz<-{cfg.comp_rz_max_abs:.2f} ({row['completeness_rz']:.2f})")
    if pd.notna(row.get("contamination_rz")) and row["contamination_rz"] > cfg.cont_rz_max:
        reasons.append(f"contamination_rz>{cfg.cont_rz_max:.2f} ({row['contamination_rz']:.2f})")

    return (len(reasons) > 0, reasons)


def act_on_folder(folder: str, apply: bool, delete: bool, quarantine_dir: Optional[str]) -> str:
    """Move to quarantine or delete; return action string."""
    if not apply:
        return "DRY-RUN"
    if delete:
        shutil.rmtree(folder)
        return "DELETED"
    else:
        assert quarantine_dir is not None
        os.makedirs(quarantine_dir, exist_ok=True)
        dest = os.path.join(quarantine_dir, os.path.basename(folder))
        # If already exists, append suffix
        i = 1
        final = dest
        while os.path.exists(final):
            final = f"{dest}__dup{i}"
            i += 1
        shutil.move(folder, final)
        return f"MOVED->{final}"


def main():
    p = argparse.ArgumentParser(
        description="Filter genome folders using robust stats + biological guards.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--qc", required=True, help="QC table CSV/TSV path.")
    p.add_argument("--base", required=True, help="Base dir containing genome subfolders.")
    p.add_argument("--out-table", default=None, help="Write annotated QC table with scores here.")
    p.add_argument("--sep", default=None, help="Field sep; auto-detect by extension if omitted.")
    p.add_argument("--apply", action="store_true", help="Actually move/delete failing genomes.")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--quarantine", default="quarantine_low_quality", help="Where to move failing genomes (default).")
    g.add_argument("--delete", action="store_true", help="Permanently delete failing genomes (DANGER).")
    # Hard-guard thresholds
    p.add_argument("--comp-min", type=float, default=90.0, help="Minimum completeness.")
    p.add_argument("--cont-max", type=float, default=5.0, help="Maximum contamination.")
    # Robust z thresholds
    p.add_argument("--num-contigs-rz-max", type=float, default=2.5, help="Max robust z for num_contigs.")
    p.add_argument("--N50-rz-max-abs", type=float, default=2.5, help="Max |rz| for low N50 (uses negative tail).")
    p.add_argument("--gc-rz-max-abs", type=float, default=3.5, help="Max |rz| for GC outliers.")
    p.add_argument("--comp-rz-max-abs", type=float, default=2.5, help="Max |rz| for low completeness (negative tail).")
    p.add_argument("--cont-rz-max", type=float, default=2.5, help="Max rz for high contamination (positive tail).")
    args = p.parse_args()

    # Auto separator
    if args.sep is None:
        if args.qc.lower().endswith(".tsv"):
            sep = "\t"
        else:
            sep = ","
    else:
        sep = args.sep

    df = pd.read_csv(args.qc, sep=sep)
    must_have_columns(df, [
        "genome", "num_contigs", "N50", "gc_percent", "completeness", "contamination",
        "fasta_path", "gbff_path"
    ])

    # Compute robust z-scores
    df["num_contigs_rz"]   = robust_z(df["num_contigs"])
    df["N50_rz"]           = robust_z(df["N50"])
    df["gc_percent_rz"]    = robust_z(df["gc_percent"])
    df["completeness_rz"]  = robust_z(df["completeness"])
    df["contamination_rz"] = robust_z(df["contamination"])

    # Evaluate pass/fail
    class Cfg: pass
    cfg = Cfg()
    cfg.comp_min = args.comp_min
    cfg.cont_max = args.cont_max
    cfg.num_contigs_rz_max = args.num_contigs_rz_max
    cfg.N50_rz_max_abs = args.N50_rz_max_abs
    cfg.gc_rz_max_abs = args.gc_rz_max_abs
    cfg.comp_rz_max_abs = args.comp_rz_max_abs
    cfg.cont_rz_max = args.cont_rz_max

    fails = []
    reasons_all = []
    for i, row in df.iterrows():
        failed, reasons = explain_fail(row, cfg)
        fails.append(failed)
        reasons_all.append("; ".join(reasons))
    df["fail"] = fails
    df["fail_reasons"] = reasons_all

    # Resolve folders and act
    base_dir = os.path.realpath(args.base)
    if not os.path.isdir(base_dir):
        raise SystemExit(f"[ERROR] Base dir does not exist: {base_dir}")

    if args.delete:
        quarantine_dir = None
    else:
        quarantine_dir = (
            os.path.realpath(args.quarantine) if args.quarantine else
            os.path.join(base_dir, "quarantine_low_quality")
        )

    actions = []
    missing_dirs = []

    for i, row in df[df["fail"]].iterrows():
        folder = safe_genome_folder(base_dir, row["genome"], row.get("gbff_path"))
        if folder is None or (not os.path.isdir(folder)):
            missing_dirs.append(row["genome"])
            actions.append("MISSING_DIR")
            continue

        # Safety: folder must be under base_dir
        real_folder = os.path.realpath(folder)
        common = os.path.commonpath([real_folder, base_dir])
        if common != base_dir:
            actions.append("SKIPPED_UNSAFE_PATH")
            continue

        action = act_on_folder(real_folder, apply=args.apply, delete=args.delete,
                               quarantine_dir=quarantine_dir)
        actions.append(action)

    df.loc[df["fail"], "action"] = actions + [""] * (df["fail"].sum() - len(actions))
    df["action"] = df.get("action", "").fillna("")

    # Summary
    total = len(df)
    n_fail = int(df["fail"].sum())
    n_pass = total - n_fail
    print("\n=== QUALITY FILTER SUMMARY ===")
    print(f"Total genomes: {total}")
    print(f"Pass: {n_pass}")
    print(f"Fail: {n_fail}")
    if missing_dirs:
        print(f"Note: {len(missing_dirs)} failing genomes had no resolvable folder under base_dir.")
    if not args.apply:
        print("\nDRY RUN (no changes made). Use --apply to move/delete failing folders.")
    else:
        if args.delete:
            print("\nApplied: failing folders were DELETED.")
        else:
            print(f"\nApplied: failing folders were MOVED to: {quarantine_dir}")

    # Save annotated table if requested (very helpful!)
    if args.out_table:
        out_sep = "\t" if args.out_table.lower().endswith(".tsv") else ","
        df.to_csv(args.out_table, index=False, sep=out_sep)
        print(f"[OK] Wrote annotated table: {args.out_table}")


if __name__ == "__main__":
    main()
