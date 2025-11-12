# E. coli Pangenome-Scale Metabolic Model Reconstruction and Analysis

**“Pangenome-Scale Metabolic Network Reconstruction Reveals a Diverse Genetic Basis for Metabolic Reactions.”**

This repository contains data and code for the pangenome-scale metabolic reconstruction of *Escherichia coli*, enabling detailed analysis of its metabolic diversity. By analyzing gene–protein–reaction (GPR) associations across thousands of strains, this work quantifies the genetic basis of *E. coli* metabolism and the evolutionary dynamics shaping metabolic functions.

> **Data availability (Zenodo)**

> * Core bundle (scripts + supporting data + panGEM): https://zenodo.org/records/17581962


---

## Method Overview (Figures)

<table>
  <tr>
    <td width="50%">
      <img src="docs/workflow.png" alt="Automated pipeline overview" width="100%">
      <p align="center"><em>Rare Metabolic Gene Essentiality is a Determinant of Microniche Adaptation in
Escherichia coli.</em></p>
    </td>
    <td width="50%">
      <img src="docs/panel1.jpg" alt="Pangenome annotation and gene-based analyses" width="100%">
      <p align="center"><em>Annotating the Pangenome Reveals the Diversity in the Genetic Basis for Metabolic Enzymes.</em></p>
    </td>
  </tr>
</table>

---

## Overview

We constructed strain-specific genome-scale metabolic models (GEMs) for **2,377 complete *E. coli* genomes**, covering ~2,700 reactions. The resulting pangenome-scale model (“panGEM”) links genotype to phenotype, supports gap-filling and targeted curation, and enables downstream analyses such as Biolog growth prediction, gene-neighborhood mining, and GPR/allele mapping.

---

## Quick Start

### 1) Create environment

```bash
# Clone
git clone https://github.com/omidard/EcopanGEM.git
cd EcopanGEM

# Conda env
conda env create -f environment.yml
conda activate panGEM

# External binaries (required)
conda install -c bioconda blast cd-hit   # provides makeblastdb/blastp/blastn + cd-hit
```

### 2) Place data (from Zenodo)

* Core bundle + GEMs: [https://zenodo.org/records/17581962)  ([Zenodo][1])

Extract into `data/` and keep this layout:

```
data/
├─ ref_model_dir/marlbr2.mat
├─ target_genome_dir/*.gbk
├─ prokka_genomes/*.gbk
├─ biolog_panGEM.csv
├─ iB21_1397.json
├─ merged_df.csv
├─ neighbors/
│  ├─ pa_gene_allt_copy2.pkl
│  ├─ locustags_genes_mapping.pkl
│  └─ pangenome_df.pkl
└─ (other folders are created by scripts)
```

> If your local paths differ from the defaults, edit the corresponding script arguments or variables in `scripts/`.

#### 2b) (Optional one-liner) Fetch core data via `Makefile`

If you prefer, you can download the core Zenodo bundle into `data/` with a single command (creates `data/Data.zip`, ignored by Git). The archive corresponds to the “Data.zip” file listed on the Zenodo record. ([Zenodo][1])

```bash
make data      # downloads to data/Data.zip
make clean     # removes downloaded data
```

---

## End-to-End Pipeline (required order)

> **Important:** The QC/QA + pangenome steps must be run **before** notebooks and downstream analyses.

### A) Genome QC & Filtering

1. **QC & (optionally) CheckM2**

   ```bash
   python scripts/genome_qc.py \
     --base_dir data/target_genome_dir \
     --out data/qc/genome_qc_results.csv \
     --threads 32 \
     --runner micromamba \
     --env ckm2 \
     --checkm2_db /path/to/checkm2_db    # omit to skip CheckM2
   ```

2. **Filter genomes (dry-run → apply)**

   ```bash
   # dry-run (prints reasons; no changes)
   python scripts/genome_qc_filtering.py \
     --qc data/qc/genome_qc_results.csv \
     --base data/target_genome_dir \
     --out-table data/qc/genome_qc_results_annotated.csv

   # apply (moves failing genomes to quarantine)
   python scripts/genome_qc_filtering.py \
     --qc data/qc/genome_qc_results.csv \
     --base data/target_genome_dir \
     --apply \
     --quarantine data/target_genome_dir/quarantine_low_quality
   ```

*(A convenience wrapper exists at `scripts/genome_qc.sh`.)*

---

### B) Pangenome Generation (GBK-aware; multi-threshold CD-HIT)

```bash
python scripts/pangenome_sensitivity_gbk.py
# Interactive prompts will ask for:
# • BASE dir (for logs/outputs), ANNOTATIONS dir, PANGENOME dir
# • cd-hit path (e.g., /usr/bin/cd-hit), threads (e.g., 32)
# • thresholds: 65,70,75,80,85,90,95  (Enter for 80)
# • -aL coverage for longer sequence (e.g., 80)
```

**Outputs (per X):**
`pangenome/cdhit_full_X.faa`, `pangenome/cdhit_full_X.faa.clstr`,
`pangenome/presence_absence_matrix_X.csv`, `pangenome/cluster_to_locus_X.json`.

**Sensitivity figures/tables**

```bash
python scripts/pangenome_sensitivity_results.py
# reads the above outputs and produces:
# figures/clusters_vs_cdhit_threshold.(png|svg)
# figures/CAR_*_cdhit80.(png|svg)  + CSV tables
```

---

### C) GEM Drafting → Gap-filling → Curation → Phenotype

1. **Orthology → draft GEMs**

   ```bash
   python scripts/pangem.py
   ```

2. **Essential gap-filling (M9-like)**

   ```bash
   python scripts/ecoli_gapfilling6.py
   ```

3. **Curation sweep (IDs, bounds, GPR merges, add/remove reactions)**

   ```bash
   python scripts/Eco_panGEM_curation.py
   ```

4. **Add spontaneous reactions (`s0001`)**

   ```bash
   python scripts/add_spont_to_gems_ecoli4.py
   ```

5. **Biolog phenotype prediction (M9 + carbon sources)**

   ```bash
   python scripts/biolog_ecoli_prediction.py
   ```

6. **Missing genes → KEGG/EC → add curated reactions (parallel)**

   ```bash
   # 4-step track:
   python scripts/gems_missing_genes.py
   python scripts/miss_locci_df_preparation.py
   python scripts/gems_missing_reactions_kegg.py
   python scripts/add_missed_reactions_ecoli.py
   ```

7. **Allele/GPR export (reaction × GEM)**

   ```bash
   python scripts/eco_gems_allels.py
   ```

8. **Gene-neighborhood mining (GBK parsing; multiprocessing)**

   ```bash
   # prepare & expand
   python scripts/genes_neighborhood_analysis_total_preparation.py
   # legacy/alternate expansion pipeline
   python scripts/genes_neighborhood_analysis_total2.py
   ```

---

## Outputs (where to find things)

* Draft models: `data/output_models_dir/`
* Gapfilled models: `data/gapfilled3/`
* Curated models: `data/gapfilled_curated/`
* Orthology/GPR tables: `data/present_absence_dir/`
* Biolog predictions: `data/biolog_data_with_predictions_panGEM_paper.csv`
* GPR/allele matrix: `data/ecoli_gprs.csv`
* Neighborhood results: `data/neighbors/`

---

## Scripts (what each does)

* `genome_qc.py` — GBFF/GenBank-aware QC; writes FASTA, contig/N50/GC stats; optional **CheckM2** completeness/contamination.
* `genome_qc_filtering.py` — robust filtering (hard guards + median/MAD) with dry-run/apply (quarantine or delete).
* `genome_qc.sh` — wrapper to run QC + CheckM2 with consistent args.
* `pangenome_sensitivity_gbk.py` — merges proteins (prefers GBK), runs **CD-HIT** at one/many identities, writes `cdhit_full_*.faa`, `*.clstr`, presence/absence matrices, and cluster→locus maps.
* `pangenome_sensitivity_results.py` — generates publication-ready figures and CSVs for cluster counts vs threshold and CAR sensitivity (@80%).
* `pangem.py` — reciprocal BLAST + nucleotide checks → presence/absence → **draft GEMs**.
* `ecoli_gapfilling6.py` — adds only **essential** missing reactions from reference under M9 conditions.
* `Eco_panGEM_curation.py` — targeted curation: ID changes, bounds, GPR merges, add/remove reactions per instruction map.
* `add_spont_to_gems_ecoli4.py` — adds spontaneous (`s0001`) reactions.
* `biolog_ecoli_prediction.py` — FBA growth predictions on M9 with specified carbon sources; appends to Biolog table.
* `gems_missing_genes.py → miss_locci_df_preparation.py → gems_missing_reactions_kegg.py → add_missed_reactions_ecoli.py` — detects genes present in genomes but absent from GEMs; maps to KEGG/EC; fetches reaction metadata/PMIDs; adds curated reactions/GPRs in parallel.
* `eco_gems_allels.py` — exports reaction×GEM allele/GPR matrix.
* `genes_neighborhood_analysis_total_preparation.py` — builds reaction-specific gene-neighborhood tables (multiprocessing; GBK parsing).
* `genes_neighborhood_analysis_total2.py` — alternate/legacy neighborhood expansion/extraction.

> Script list cross-checked with your repo’s `scripts/` directory. ([GitHub][3])

---

## Environment & Requirements

* **Python:** ≥ 3.10
* **Packages:** see `environment.yml` in this repo (COBRApy, pandas, numpy, scipy, tqdm, Biopython, bioservices, matplotlib, rich, requests, joblib, scikit-learn, openpyxl, pyyaml, python-dotenv, etc.).
* **External binaries:** **BLAST+** (`makeblastdb`, `blastp`, `blastn`) and **CD-HIT** must be on `PATH`.
* **Hardware:** Several steps default to high parallelism (64–96 workers). Reduce pool sizes if RAM/CPU is constrained.

---

## Notes & Tips

* Historical absolute paths (e.g., `/home/omidard/...`) are mirrored by the `data/` layout above. If you use a different structure, pass explicit CLI args or edit the variables at the top of the scripts.
* KEGG queries via `bioservices` may rate-limit; consider chunking or local caching.
* Some curation code expects `.json` model files and specific IDs/bounds (kept consistent across steps here).

---

## Citation

If you use this repository or derived models, please cite:

> **“Annotating the Pangenome Reveals the Diversity in the Genetic Basis for Metabolic Enzymes”**

---

## Contact

**Omid Ardalani** — [omidard@biosustain.dtu.dk](mailto:omidard@biosustain.dtu.dk) · Issues and PRs welcome.

---
