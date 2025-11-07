#!/usr/bin/env bash
# panGEM: end-to-end runner (copy this into README as a single block)
# Purpose
# - One command to execute the E. coli panGEM pipeline with sensible checks.
# - Runs each script in the expected order and prints where outputs land.
#
# Usage
#   bash run_pangem.sh
#
# Requirements
#   - Python 3.10+ with: cobra, pandas, biopython (Bio), tqdm, bioservices, openpyxl
#   - NCBI BLAST+ on PATH: makeblastdb, blastp, blastn
#   - Data from Zenodo placed under ./data (see paths below)

set -euo pipefail

# --------- CONFIG (edit as needed) ---------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${REPO_ROOT}/scripts"
DATA="${REPO_ROOT}/data"

# Core inputs
REF_MODEL="${DATA}/ref_model_dir/marlbr2.mat"
TARGET_GBKS_DIR="${DATA}/target_genome_dir"     # *.gbk for strains
PROKKA_GBKS_DIR="${DATA}/prokka_genomes"       # *.gbk (if used by neighborhood step)
DONOR_MODEL="${DATA}/iB21_1397.json"           # optional, used by Eco_panGEM_curation.py
BIOLOG_CSV="${DATA}/biolog_panGEM.csv"

# Intermediate/output dirs (created by scripts)
OUTPUT_MODELS_DIR="${DATA}/output_models_dir"
GAPFILLED_DIR="${DATA}/gapfilled"
GAPFILLED3_DIR="${DATA}/gapfilled3"
CURATED_DIR="${DATA}/gapfilled_curated"
NEI_DIR="${DATA}/neighbors"

# Optional inputs for neighborhood expansion
PA_PICKLE="${NEI_DIR}/pa_gene_allt_copy2.pkl"
MAP_PICKLE="${NEI_DIR}/locustags_genes_mapping.pkl"
PG_PKL="${NEI_DIR}/pangenome_df.pkl"

# Optional inputs for allele mapping
ALLELE_MAP="${DATA}/header_to_allele.pickle"
COMPLETE_GEMS_CSV="${REPO_ROOT}/complete_gems.csv"

# Optional curated-reaction table step
MERGED_TABLE="${DATA}/merged_df.csv"   # used by add_missed_reactions_ecoli.py

# Python interpreter
PYTHON_BIN="${PYTHON_BIN:-python}"

# ------------------------------------------

say() { printf "\n[%s] %s\n" "$(date +%H:%M:%S)" "$*"; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1"; exit 1; }; }

# 0) Dependency & data checks
say "Checking dependencies"
need "${PYTHON_BIN}"
need makeblastdb
need blastp
need blastn

"${PYTHON_BIN}" - <<'PY'
import sys
mods = ["cobra","pandas","Bio","tqdm","openpyxl","bioservices"]
missing = []
for m in mods:
    try: __import__(m)
    except Exception: missing.append(m)
if missing:
    print("Missing Python modules:", ", ".join(missing))
    sys.exit(1)
PY

say "Checking data"
test -f "${REF_MODEL}" || { echo "Not found: ${REF_MODEL}"; exit 1; }
test -d "${TARGET_GBKS_DIR}" || { echo "Not found: ${TARGET_GBKS_DIR}"; exit 1; }
ls -1 "${TARGET_GBKS_DIR}"/*.gbk >/dev/null 2>&1 || { echo "No .gbk in ${TARGET_GBKS_DIR}"; exit 1; }

mkdir -p "${GAPFILLED_DIR}" "${GAPFILLED3_DIR}" "${CURATED_DIR}" "${NEI_DIR}"

# 1) Orthology & draft GEMs
if [ -f "${SCRIPTS}/pangem.py" ]; then
  say "1/10 pangem.py → draft GEMs & orthology matrices"
  ( cd "${DATA}" && "${PYTHON_BIN}" "${SCRIPTS}/pangem.py" )
  say "   outputs: present_absence_dir/*, initial_models_dir/*.json, ${OUTPUT_MODELS_DIR}/*.json"
else
  say "1/10 pangem.py skipped (not found: ${SCRIPTS}/pangem.py)"
fi

# 2) Genome CDS vs GEM genes (missing loci)
if [ -f "${SCRIPTS}/gems_missing_genes.py" ]; then
  say "2/10 gems_missing_genes.py → ecoli_missing_locus_tags.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/gems_missing_genes.py"
else
  say "2/10 gems_missing_genes.py skipped"
fi

# 3) Tidy missing loci
if [ -f "${SCRIPTS}/miss_locci_df_preparation.py" ]; then
  say "3/10 miss_locci_df_preparation.py → tidy_missing_locci_ecoli.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/miss_locci_df_preparation.py"
else
  say "3/10 miss_locci_df_preparation.py skipped"
fi

# 4) KEGG mapping + EC metadata
if [ -f "${SCRIPTS}/gems_missing_reactions_kegg.py" ]; then
  say "4/10 gems_missing_reactions_kegg.py → ecoli_missed_reactions_refined.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/gems_missing_reactions_kegg.py"
else
  say "4/10 gems_missing_reactions_kegg.py skipped"
fi

# 5) Add/modify reactions from curated table (optional; requires MERGED_TABLE)
if [ -f "${SCRIPTS}/add_missed_reactions_ecoli.py" ] && [ -f "${MERGED_TABLE}" ]; then
  say "5/10 add_missed_reactions_ecoli.py → gapfilled3/*.json (from ${MERGED_TABLE})"
  "${PYTHON_BIN}" "${SCRIPTS}/add_missed_reactions_ecoli.py"
else
  say "5/10 add_missed_reactions_ecoli.py skipped (script or ${MERGED_TABLE} not found)"
fi

# 6) Gap-filling under M9 (essential subset)
if [ -f "${SCRIPTS}/ecoli_gapfilling6.py" ]; then
  say "6/10 ecoli_gapfilling6.py → gapfilled3/*.json"
  "${PYTHON_BIN}" "${SCRIPTS}/ecoli_gapfilling6.py"
else
  say "6/10 ecoli_gapfilling6.py skipped"
fi

# 7) Add “spontaneous” (s0001) reactions
if [ -f "${SCRIPTS}/add_spont_to_gems_ecoli4.py" ]; then
  say "7/10 add_spont_to_gems_ecoli4.py → gapfilled/*.json (in-place)"
  "${PYTHON_BIN}" "${SCRIPTS}/add_spont_to_gems_ecoli4.py"
else
  say "7/10 add_spont_to_gems_ecoli4.py skipped"
fi

# 8) Targeted curation pass
if [ -f "${SCRIPTS}/Eco_panGEM_curation.py" ]; then
  if [ -f "${DONOR_MODEL}" ]; then
    say "8/10 Eco_panGEM_curation.py → gapfilled_curated/*.json.json"
    "${PYTHON_BIN}" "${SCRIPTS}/Eco_panGEM_curation.py"
  else
    say "8/10 Eco_panGEM_curation.py skipped (donor model missing: ${DONOR_MODEL})"
  fi
else
  say "8/10 Eco_panGEM_curation.py skipped"
fi

# 9) Biolog phenotype prediction on M9
if [ -f "${SCRIPTS}/biolog_ecoli_prediction.py" ] && [ -f "${BIOLOG_CSV}" ]; then
  say "9/10 biolog_ecoli_prediction.py → biolog_data_with_predictions_panGEM_paper.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/biolog_ecoli_prediction.py"
else
  say "9/10 biolog_ecoli_prediction.py skipped (script or ${BIOLOG_CSV} missing)"
fi

# 10) Reaction×GEM GPR matrix (alleles)
if [ -f "${SCRIPTS}/eco_gems_allels.py" ] && [ -f "${ALLELE_MAP}" ] && [ -f "${COMPLETE_GEMS_CSV}" ]; then
  say "10/10 eco_gems_allels.py → ecoli_gprs.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/eco_gems_allels.py"
else
  say "10/10 eco_gems_allels.py skipped (script or inputs missing)"
fi

# Neighborhood analysis — Stage A (expand P/A → rows)
if [ -f "${SCRIPTS}/genes_neighborhood_analysis_total_preparation.py" ] && \
   [ -f "${PA_PICKLE}" ] && [ -f "${MAP_PICKLE}" ] && [ -f "${PG_PKL}" ]; then
  say "Neighborhood A: expand presence/absence → final_gene_neighborhood.csv"
  "${PYTHON_BIN}" "${SCRIPTS}/genes_neighborhood_analysis_total_preparation.py"
else
  say "Neighborhood A skipped (script or pickles missing)"
fi

# Neighborhood analysis — Stage B (extract neighbors from GenBank)
# NOTE: If you keep both “expand” and “extract” scripts under the same name, only one will run.
# Rename the extractor to genes_neighborhood_extract.py for clarity and enable below.
EXTRACTOR="${SCRIPTS}/genes_neighborhood_extract.py"
if [ -f "${EXTRACTOR}" ] && [ -d "${PROKKA_GBKS_DIR}" ]; then
  say "Neighborhood B: extract flanking CDS annotations → all_reactions_gene_neighborhood.csv"
  "${PYTHON_BIN}" "${EXTRACTOR}"
else
  say "Neighborhood B skipped (genes_neighborhood_extract.py or GenBank dir missing)"
fi

# Summary
say "DONE."
echo " Key locations:"
echo "  - Draft models:        ${OUTPUT_MODELS_DIR}"
echo "  - Gapfilled models:    ${GAPFILLED3_DIR}"
echo "  - Curated models:      ${CURATED_DIR}"
echo "  - Orthology matrices:  ${DATA}/present_absence_dir"
echo "  - Biolog predictions:  ${DATA}/biolog_data_with_predictions_panGEM_paper.csv (if run)"
echo "  - GPR matrix:          ${DATA}/ecoli_gprs.csv (if run)"
echo "  - Neighborhood tables: ${NEI_DIR} (if run)"
