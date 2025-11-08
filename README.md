Here’s a cleaned, drop-in replacement of your current README content. I focused on consistency, correctness, and reviewer-friendly clarity (kept your structure, trimmed ambiguity, fixed thresholds/typos, and standardized terminology and paths).

---

# E. coli Pangenome-Scale Metabolic Model Reconstruction and Analysis: “Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions”

This repository contains data and code for the pangenome-scale metabolic reconstruction of *Escherichia coli*, enabling detailed analysis of its metabolic diversity. By analyzing gene–protein–reaction (GPR) associations across thousands of strains, this work quantifies the genetic basis of *E. coli* metabolism and the evolutionary dynamics shaping metabolic functions.

## Method Overview (Figures)

<table>
  <tr>
    <td width="50%">
      <img src="docs/workflow.png" alt="Automated pipeline overview" width="100%">
      <p align="center"><em>Pipeline overview for panGEM reconstruction and analysis.</em></p>
    </td>
    <td width="50%">
      <img src="docs/panel1.jpg" alt="Pangenome annotation and gene-based analyses" width="100%">
      <p align="center"><em>Pangenome annotation and gene-based analyses.</em></p>
    </td>
  </tr>
</table>

## Overview

We constructed strain-specific genome-scale metabolic models (GEMs) for 2,377 complete *E. coli* genomes, covering ~2,700 reactions. The resulting pangenome-scale model (“panGEM”) links genotype to phenotype, supports gap-filling and targeted curation, and enables downstream analyses such as Biolog growth prediction, gene-neighborhood mining, and GPR/allele mapping.

## Quick Start

1. **Create environment**

   ```bash
   conda env create -f environment.yml
   conda activate panGEM
   conda install -c bioconda blast cd-hit   # makeblastdb/blastp/blastn + cd-hit
   ```
2. **Place data (from Zenodo)**

   * Core bundle: [https://zenodo.org/records/14028473](https://zenodo.org/records/14028473)
   * GEMs bundle: [https://zenodo.org/records/13825392](https://zenodo.org/records/13825392)

   Extract into `data/` keeping this layout:

   ```
   data/
   ├─ ref_model_dir/marlbr2.mat
   ├─ target_genome_dir/*.gbk
   ├─ prokka_genomes/*.gbk
   ├─ biolog_panGEM.csv
   ├─ iB21_1397.json
   ├─ merged_df.csv                 # if using curated-reaction addition
   ├─ neighbors/
   │  ├─ pa_gene_allt_copy2.pkl
   │  ├─ locustags_genes_mapping.pkl
   │  └─ pangenome_df.pkl
   └─ (other folders are created by scripts)
   ```
3. **Run core pipeline**

   ```bash
   python scripts/pangem.py                    # orthology → draft GEMs
   python scripts/ecoli_gapfilling6.py         # essential gap-filling (M9-like)
   python scripts/Eco_panGEM_curation.py       # targeted curation sweep
   python scripts/biolog_ecoli_prediction.py   # phenotype (Biolog/M9) prediction
   ```

   *Optional steps are listed below; script paths assume the layout above.*

## Repository Structure

* **data/** — input data and generated outputs (see layout above).
* **scripts/** — standalone steps for **QC → pangenome → panGEM → analyses**

  * `genome_qc.py` — GBFF/GenBank–aware genome QC; exports FASTA, assembly stats, and optional **CheckM2** completeness/contamination.
  * `genome_qc_filtering.py` — robust filtering of low-quality genomes (hard guards + median/MAD); dry-run by default; can quarantine or delete.
  * `genome_qc.sh` — convenience wrapper to run QC + CheckM2 with consistent args.
  * `pangenome_sensitivity_gbk.py` — GBK-aware protein merge and **CD-HIT** clustering at one or multiple identity thresholds; writes `cdhit_full_*.faa`, `*.clstr`, `presence_absence_matrix_*.csv`, `cluster_to_locus_*.json`.
  * `pangenome_sensitivity_results.py` — figures/CSV tables for cluster counts vs threshold and CAR (Core/Accessory/Rare) sensitivity at 80%.
  * `pangem.py` — reciprocal BLAST + nucleotide checks → presence/absence matrices; drafts strain-specific GEMs.
  * `ecoli_gapfilling6.py` — adds essential subset from reference under M9 conditions (keeps only missing reactions essential for growth).
  * `Eco_panGEM_curation.py` — targeted curation (ID changes, bounds, GPR merges, add/remove reactions per instruction map).
  * `add_spont_to_gems_ecoli4.py` — adds spontaneous (`s0001`) reactions to models (optional).
  * `biolog_ecoli_prediction.py` — FBA growth predictions on M9 with carbon sources; writes predictions back to the Biolog table.
  * `gems_missing_genes.py` → `miss_locci_df_preparation.py` → `gems_missing_reactions_kegg.py` → `add_missed_reactions_ecoli.py` — optional track to:

    * detect genes present in genomes but absent from GEMs,
    * parse/clean product strings; separate hypotheticals,
    * map candidates to KEGG/EC and fetch reaction metadata/PMIDs,
    * add curated reactions/GPRs back to models in parallel.
  * `eco_gems_allels.py` — reaction×GEM matrix of allele/GPR content (export).
  * `genes_neighborhood_analysis_total_preparation.py` — builds reaction-specific gene-neighborhood tables across genomes (multiprocessing; GBK parsing).
  * `genes_neighborhood_analysis_total2.py` — alternative/legacy implementation of neighborhood expansion/extraction.

## Typical Outputs

* Draft models: `data/output_models_dir/`
* Gapfilled models: `data/gapfilled3/`
* Curated models: `data/gapfilled_curated/`
* Orthology/GPR tables: `data/present_absence_dir/`
* Biolog predictions: `data/biolog_data_with_predictions_panGEM_paper.csv`
* GPR/allele matrix: `data/ecoli_gprs.csv`
* Neighborhood results: `data/neighbors/`

## Requirements

* **Software:** Python ≥ 3.10, COBRApy, pandas, Biopython, tqdm, openpyxl, bioservices; BLAST+ and CD-HIT on `PATH`.
* **Hardware:** Several steps default to high parallelism (64–96 workers). Reduce pool sizes if RAM/CPU is constrained.

## Notes & Tips

* Some scripts historically used absolute paths (e.g., `/home/omidard/...`). The `data/` layout mirrors those expectations; adjust paths inside scripts if you use a different structure.
* KEGG queries via `bioservices` can rate-limit; consider caching.
* File names ending with `.json` are expected by `Eco_panGEM_curation.py`.

---

## Data (downloaded from Zenodo on demand)

This repo keeps large data **out of Git history**. A small `Makefile` downloads the archive when you need it.

**What you get**

* Source: Zenodo record `14028473`
* File: `data/Data.zip` (created locally on demand; ignored by Git)

**Requirements**

* `make`, `curl`

**Usage**

```bash
make data
```

The file appears at `data/Data.zip`.

**Clean up**

```bash
make clean
```

**Notes**

* `data/` is in `.gitignore`, so you won’t accidentally commit large files.
* On Windows, use Git Bash or WSL to run `make` and `curl`.

## Genome QC/QA and Pangenome Generation (run **before** the notebooks)

Minimal, reproducible steps to (i) QC/QA genomes, (ii) filter low-quality assemblies, and (iii) build the pangenome used downstream.

### 0) Prerequisites

* **Inputs:** per-genome annotation folders each with a single GenBank file (`*.gbff`/`*.gbk`) *or* a flat folder of `*.gbk`.
* **Tools:** `python>=3.10`, `biopython`, `pandas`, `tqdm`, `matplotlib`, `cd-hit`; optional **CheckM2** via `micromamba/conda`.

  ```bash
  micromamba create -n ckm2 -c conda-forge -c bioconda checkm2 biopython pandas tqdm matplotlib
  ```

### 1) QC each genome (GBFF-aware) and run CheckM2

**Script:** `scripts/genome_qc.py`
**Input:** `--base_dir` with subfolders per genome, each containing one `*.gbff/*.gbk`
**Output:** per-genome QC table with contig stats + optional completeness/contamination

```bash
python scripts/genome_qc.py \
  --base_dir /path/to/annotations_root \
  --out /path/to/qc/genome_qc_results.csv \
  --threads 32 \
  --runner micromamba \
  --env ckm2 \
  --checkm2_db /path/to/checkm2_db
```

*(Omit `--checkm2_db` to skip CheckM2; columns will be NA.)*

### 2) Filter genomes by quality (dry-run → apply)

**Script:** `scripts/genome_qc_filtering.py`
**Input:** QC table from step 1 and the **same base directory** of genome folders
**Action:** mark failing genomes via hard guards and robust Z-scores; dry-run by default

```bash
# Dry run (inspect summary and reasons per genome)
python scripts/genome_qc_filtering.py \
  --qc /path/to/qc/genome_qc_results.csv \
  --base /path/to/annotations_root \
  --out-table /path/to/qc/genome_qc_results_annotated.csv

# Apply: move failing folders to quarantine (recommended)
python scripts/genome_qc_filtering.py \
  --qc /path/to/qc/genome_qc_results.csv \
  --base /path/to/annotations_root \
  --apply \
  --quarantine /path/to/annotations_root/quarantine_low_quality
```

**Defaults:** `completeness ≥ 90`, `contamination ≤ 5` (tune with `--comp-min/--cont-max`).
**Danger:** `--delete` permanently removes failing folders (mutually exclusive with `--quarantine`).

### 3) Build the pangenome (CD-HIT; GBK-aware; multi-threshold)

**Script:** `scripts/pangenome_sensitivity_gbk.py`
**Input:** an **annotations** directory containing either a flat collection of `*.gbk` or a per-genome tree with `*.faa`.
**Outputs (per threshold X):**
`pangenome/cdhit_full_X.faa`, `pangenome/cdhit_full_X.faa.clstr`,
`pangenome/presence_absence_matrix_X.csv`, `pangenome/cluster_to_locus_X.json`.

```bash
python scripts/pangenome_sensitivity_gbk.py
# Prompts:
#  • BASE dir (logs/outputs), ANNOTATIONS dir, PANGENOME dir
#  • path to cd-hit (e.g., /usr/bin/cd-hit)
#  • threads (e.g., 32)
#  • identity thresholds, e.g.: 65,70,75,80,85,90,95   (Enter for 80)
#  • alignment coverage for longer sequence (-aL), e.g., 80
```

**What later notebooks expect (defaults):** `presence_absence_matrix_80.csv`, `cluster_to_locus_80.json` under the chosen `pangenome/` directory.

### 4) Summarize pangenome sensitivity (figures + tables)

**Script:** `scripts/pangenome_sensitivity_results.py`
**Input:** the JSON/CSV outputs from step 3 under a single `BASE` directory
**Outputs:** publication-ready panels (`.png/.svg`) and CSV tables in `BASE/figures/`

```bash
python scripts/pangenome_sensitivity_results.py
```

## Citation

If you use this repository or derived models, please cite:
**“Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions.”**
Zenodo datasets: **14028473** (core) and **13825392** (GEMs).

## Contact

Omid Ardalani — [omidard@biosustain.dtu.dk](mailto:omidard@biosustain.dtu.dk) · Issues and PRs welcome.
