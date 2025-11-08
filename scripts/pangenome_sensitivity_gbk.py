#!/usr/bin/env python3
"""
pangenome_analysis_minimal.py  —  GBK-aware, multi-threshold CD-HIT

WHAT CHANGED (addresses your points):
1) Multi-threshold input:
   - You can now enter: 65,70,75,80,85,90,95   OR   [65, 70, 75, 80, 85, 90, 95]
   - The script will loop over each threshold and produce separate outputs per threshold.
   - Why it previously only ran 80: the old code tried to `int()` your whole string
     "[65, 70, 75, 80, 85, 90, 95]" which threw a ValueError and fell back to default 80.

2) Identity vs threads (-c vs -T in cd-hit):
   - In CD-HIT, **-c** is the minimum identity (0.00–1.00). **-T** is *threads*, not identity.
   - This script sets identity with **-c**. When you enter "80", it uses `-c 0.80`.
   - The fixed alignment coverage for the longer sequence is set with **-aL** (here: 0.80).

3) Word size (-n) handling:
   - You asked for `-n 5`. CD-HIT only allows `-n 5` when identity `-c >= 0.70`.
   - For thresholds < 70 (e.g., 65), CD-HIT requires `-n 4` (and 50–60% would require `-n 3`).
   - This script **prefers `-n 5`** but automatically downgrades to the allowed word length
     when a threshold is too low for `-n 5`. A clear INFO log is emitted when it downgrades.

4) Alignment coverage (-aL):
   - `-aL` is fixed across thresholds and interpreted as a fraction (0.80 = 80%).
   - You can change the default once and it will apply to all thresholds.

INPUTS/OUTPUTS:
- ANNOTATIONS can be a flat directory of *.gbk (no subfolders), or the original per-genome *.faa tree.
- The merged proteins will always carry genome IDs in headers like:
      >{seq_id}|{genome_id}|{locus_tag}|{protein_id}
  so presence/absence per genome is unambiguous even when locus_tags don’t encode genome names.

- For each threshold X, outputs are:
    pangenome/cdhit_full_X.faa
    pangenome/cdhit_full_X.faa.clstr
    pangenome/presence_absence_matrix_X.csv
    pangenome/cluster_to_locus_X.json
"""

import os
import sys
import re
import json
import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Tuple
import pandas as pd
from Bio import SeqIO
from multiprocessing import cpu_count


# ----------------------------
# Interactive helpers
# ----------------------------
def ask_dir(prompt: str, default: Optional[Path] = None, must_exist: bool = True) -> str:
    if default is not None:
        raw = input(f"{prompt} (Enter for '{default}') : ").strip().strip('"').strip("'")
        p = Path(os.path.expanduser(raw)).resolve() if raw else default
    else:
        raw = input(prompt + " ").strip().strip('"').strip("'")
        if not raw:
            print("No path provided."); sys.exit(1)
        p = Path(os.path.expanduser(raw)).resolve()
    if must_exist and not p.is_dir():
        print(f"Not a directory: {p}"); sys.exit(1)
    if not must_exist:
        p.mkdir(parents=True, exist_ok=True)
    return str(p)

def ask_file(prompt: str, default: Optional[Path] = None) -> str:
    if default is not None:
        raw = input(f"{prompt} (Enter for '{default}') : ").strip().strip('"').strip("'")
        f = Path(os.path.expanduser(raw)).resolve() if raw else default
    else:
        raw = input(prompt + " ").strip().strip('"').strip("'")
        if not raw:
            print("No file provided."); sys.exit(1)
        f = Path(os.path.expanduser(raw)).resolve()
    return str(f)

def parse_thresholds_input(raw: str) -> List[int]:
    """
    Accepts: "65,70,75"  |  "65 70 75"  |  "[65, 70, 75]"
    Returns sorted unique integer thresholds (50..100).
    """
    s = raw.strip()
    s = s.strip("[]()")            # drop brackets if present
    s = s.replace(",", " ")        # unify separators
    parts = [p for p in s.split() if p]
    out = []
    for p in parts:
        try:
            v = int(p)
            if 50 <= v <= 100:
                out.append(v)
        except ValueError:
            pass
    out = sorted(set(out))
    return out


# ----------------------------
# Logging
# ----------------------------
def setup_logger(log_file: str) -> logging.Logger:
    logger = logging.getLogger("pangenome_minimal")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler(log_file); fh.setLevel(logging.DEBUG); fh.setFormatter(fmt)
    ch = logging.StreamHandler();      ch.setLevel(logging.INFO);  ch.setFormatter(fmt)
    if not logger.handlers:
        logger.addHandler(fh)
        logger.addHandler(ch)
    return logger


# ----------------------------
# Header / parsing helpers
# ----------------------------
def build_header(genome_id: str, locus_tag: Optional[str], protein_id: Optional[str], seq_ix: int) -> str:
    """
    Build header that *always* includes genome_id as second field:
        >{seq_id}|{genome_id}|{locus_tag}|{protein_id}
    """
    lt = locus_tag if locus_tag else "NA"
    pid = protein_id if protein_id else "NA"
    seq_id = locus_tag or protein_id or f"{genome_id}_cds{seq_ix}"
    return f"{seq_id}|{genome_id}|{lt}|{pid}"

def extract_genome_id_from_header(header: str) -> str:
    h = header.lstrip(">")
    parts = h.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    m = re.search(r'(GCA|GCF)_[0-9]+(?:\.[0-9]+)?', h)
    if m:
        return m.group(0)
    if len(parts) >= 1 and re.fullmatch(r'\d{3,}\.\d+', parts[0]):
        return parts[0]
    return "unknown"


# ----------------------------
# Merge proteins (GBK preferred)
# ----------------------------
def merge_proteins_from_faa_tree(annotation_dir: str, merged_faa: str, logger: logging.Logger) -> int:
    seqs = 0
    with open(merged_faa, 'w') as merged_file:
        for root, _, files in os.walk(annotation_dir):
            for file in files:
                if file.endswith(".faa") and "hypotheticals" not in file.lower():
                    faa_path = os.path.join(root, file)
                    genome_folder = os.path.basename(root)
                    with open(faa_path, 'r') as f:
                        for line in f:
                            if line.startswith(">"):
                                if "|" not in line:
                                    merged_file.write(line.strip() + f"|{genome_folder}|NA|NA\n")
                                else:
                                    merged_file.write(line.strip() + "\n")
                            else:
                                if line.strip():
                                    seqs += 1
                                merged_file.write(line)
    logger.info(f"Merged protein FASTA (FAA tree) saved to {merged_faa} with ~{seqs} sequences.")
    return seqs

def extract_proteins_from_gbk(gbk_path: Path, logger: logging.Logger) -> List[Tuple[str, str]]:
    genome_id = gbk_path.stem
    out = []
    cds_ix = 0
    missing = 0
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            cds_ix += 1
            quals = feat.qualifiers
            aa = quals.get("translation", [None])[0]
            if not aa:
                missing += 1
                continue
            locus_tag = quals.get("locus_tag", [None])[0]
            protein_id = quals.get("protein_id", [None])[0]
            header = build_header(genome_id, locus_tag, protein_id, cds_ix)
            out.append((header, aa.replace(" ", "").replace("\n", "")))
    if missing:
        logger.debug(f"{gbk_path.name}: CDS={cds_ix}, missing translations={missing}")
    return out

def merge_proteins_from_gbks_flat(annotation_dir: str, merged_faa: str, logger: logging.Logger) -> int:
    gbks = sorted(Path(annotation_dir).glob("*.gbk"))
    if not gbks:
        return 0
    total = 0
    with open(merged_faa, "w") as out_f:
        for i, gbk in enumerate(gbks, 1):
            try:
                recs = extract_proteins_from_gbk(gbk, logger)
                for hdr, seq in recs:
                    out_f.write(f">{hdr}\n{seq}\n")
                total += len(recs)
                if i % 50 == 0:
                    logger.info(f"Processed {i}/{len(gbks)} GBKs → {total} proteins")
            except Exception as e:
                logger.warning(f"Failed to parse {gbk}: {e}")
    logger.info(f"Merged protein FASTA (GBK flat) saved to {merged_faa} with {total} sequences.")
    return total

def merge_proteins(annotation_dir: str, merged_faa: str, logger: logging.Logger) -> None:
    if any(Path(annotation_dir).glob("*.gbk")):
        n = merge_proteins_from_gbks_flat(annotation_dir, merged_faa, logger)
        if n == 0:
            logger.error("No protein sequences extracted from GBKs."); sys.exit(1)
    else:
        n = merge_proteins_from_faa_tree(annotation_dir, merged_faa, logger)
        if n == 0:
            logger.error("No protein sequences found in FAA tree."); sys.exit(1)


# ----------------------------
# CD-HIT execution
# ----------------------------
def pick_word_length(identity_fraction: float) -> int:
    """
    CD-HIT constraints for proteins:
      c >= 0.70 -> n = 5
      0.60 <= c < 0.70 -> n = 4
      0.50 <= c < 0.60 -> n = 3
    """
    if identity_fraction >= 0.70:
        return 5
    if identity_fraction >= 0.60:
        return 4
    return 3

def run_cd_hit(
    faa_input: str,
    identity_percent: int,
    output_file: str,
    logger: logging.Logger,
    cd_hit_bin: str,
    threads: int,
    align_cov_percent: int
) -> None:
    c = max(50, min(100, identity_percent)) / 100.0
    # Prefer n=5, but honor CD-HIT constraints:
    n_allowed = pick_word_length(c)
    if n_allowed != 5:
        logger.info(f"Threshold {identity_percent}% requires -n {n_allowed}; using that instead of -n 5.")
    n = n_allowed

    aL = max(0, min(100, align_cov_percent)) / 100.0  # fraction

    cmd = [
        cd_hit_bin,
        "-i", faa_input,
        "-o", output_file,
        "-c", f"{c:.2f}",      # identity
        "-n", str(n),          # word length (constrained)
        "-aL", f"{aL:.2f}",    # alignment coverage for LONGER seq
        "-T", str(int(threads)),
        "-M", "16000",
        "-d", "0",
    ]
    logger.info(f"CD-HIT @ {identity_percent}% (c={c:.2f}, n={n}, aL={aL:.2f}) → {output_file}")
    subprocess.run(cmd, check=True)


# ----------------------------
# Parse clusters and build PA
# ----------------------------
def parse_cd_hit_clusters(cluster_file: str, logger: logging.Logger) -> dict:
    clusters = {}
    id_pattern = re.compile(r'>(\S+)')
    current = None
    with open(cluster_file) as f:
        for line in f:
            s = line.strip()
            if s.startswith(">Cluster"):
                current = s.split()[1]
                clusters[current] = []
            else:
                m = id_pattern.search(s)
                if m and current is not None:
                    hdr = m.group(1)  # e.g., seq|genome|locus|prot
                    clusters[current].append(hdr)
    logger.info(f"Parsed {len(clusters)} clusters from {cluster_file}")
    return clusters

def build_presence_absence_matrix(clusters: dict, logger: logging.Logger) -> pd.DataFrame:
    genomes = set()
    per_cluster = {}
    for cid, members in clusters.items():
        gset = set()
        for hdr in members:
            gid = extract_genome_id_from_header(hdr)
            genomes.add(gid)
            gset.add(gid)
        per_cluster[cid] = gset
    all_genomes = sorted(genomes)
    data = [[1 if g in per_cluster[cid] else 0 for g in all_genomes] for cid in per_cluster]
    df = pd.DataFrame(data, index=list(per_cluster.keys()), columns=all_genomes)
    logger.info(f"PA matrix: {df.shape[0]} clusters × {df.shape[1]} genomes")
    return df

def clusters_to_rich_mapping(clusters: dict) -> dict:
    out = {}
    for cid, members in clusters.items():
        lst = []
        for hdr in members:
            parts = hdr.lstrip(">").split("|")
            seq_id  = parts[0] if len(parts) > 0 else "NA"
            genome  = parts[1] if len(parts) > 1 else "unknown"
            locus   = parts[2] if len(parts) > 2 else "NA"
            prot_id = parts[3] if len(parts) > 3 else "NA"
            lst.append({"seq_id": seq_id, "genome_id": genome, "locus_tag": locus, "protein_id": prot_id})
        out[cid] = lst
    return out


# ----------------------------
# Main
# ----------------------------
def main():
    # Paths
    default_base = Path.cwd()
    base_dir = ask_dir("Enter BASE directory (for logs & defaults)", default=default_base, must_exist=False)
    default_annotations = Path(base_dir) / "annotations"
    default_pangenome   = Path(base_dir) / "pangenome"
    annotation_dir = ask_dir("Enter ANNOTATIONS dir (contains *.gbk OR per-genome *.faa):",
                             default=default_annotations, must_exist=True)
    pangenome_dir = ask_dir("Enter PANGENOME output dir:", default=default_pangenome, must_exist=False)
    cd_hit_bin = ask_file("Enter path to CD-HIT binary:", default=Path("/home/omidard/cdhit/cd-hit"))

    # Config
    try:
        threads_in = input(f"Threads for CD-HIT (Enter for {cpu_count()}): ").strip()
        threads = int(threads_in) if threads_in else cpu_count()
    except Exception:
        threads = cpu_count()

    # Multi-threshold parse
    thr_raw = input("CD-HIT identity threshold(s) % (e.g., 65,70,75 or [65,70,75]; Enter for 80): ").strip()
    thresholds = parse_thresholds_input(thr_raw) if thr_raw else [80]
    if not thresholds:
        thresholds = [80]

    # Fixed alignment coverage (longer sequence), percent → applies to all thresholds
    aL_raw = input("Alignment coverage for longer sequence, -aL % (Enter for 80): ").strip()
    try:
        align_cov_percent = int(aL_raw) if aL_raw else 80
    except Exception:
        align_cov_percent = 80

    # Logger
    log_file = str(Path(base_dir) / "pangenome_minimal.log")
    logger = setup_logger(log_file)
    logger.info(f"Starting pangenome pipeline (GBK-aware). Thresholds: {thresholds}, aL={align_cov_percent}%")

    # 1) Merge proteins
    merged_faa = str(Path(pangenome_dir) / "merged_proteins.faa")
    merge_proteins(annotation_dir, merged_faa, logger)
    if not Path(merged_faa).is_file() or os.stat(merged_faa).st_size == 0:
        logger.error("Merged protein FASTA file is empty or missing. Aborting.")
        sys.exit(1)

    # 2) Run CD-HIT for each threshold
    for thr in thresholds:
        out_faa = str(Path(pangenome_dir) / f"cdhit_full_{thr}.faa")
        run_cd_hit(
            faa_input=merged_faa,
            identity_percent=thr,
            output_file=out_faa,
            logger=logger,
            cd_hit_bin=cd_hit_bin,
            threads=threads,
            align_cov_percent=align_cov_percent
        )

        # 3) Parse clusters & PA for this threshold
        cluster_file = out_faa + ".clstr"
        if not Path(cluster_file).is_file():
            logger.error(f"Cluster file not found: {cluster_file}")
            continue

        clusters_raw = parse_cd_hit_clusters(cluster_file, logger)
        presence_df = build_presence_absence_matrix(clusters_raw, logger)

        # 4) Save outputs for this threshold
        presence_csv = str(Path(pangenome_dir) / f"presence_absence_matrix_{thr}.csv")
        presence_df.to_csv(presence_csv)
        logger.info(f"[{thr}%] Presence/absence matrix → {presence_csv}")

        clusters_json = str(Path(pangenome_dir) / f"cluster_to_locus_{thr}.json")
        with open(clusters_json, "w") as f:
            json.dump(clusters_to_rich_mapping(clusters_raw), f, indent=2, ensure_ascii=False)
        logger.info(f"[{thr}%] Cluster-to-locus mapping → {clusters_json}")

    logger.info("Done.")

if __name__ == "__main__":
    main()
