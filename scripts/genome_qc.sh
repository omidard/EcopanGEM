#!/usr/bin/env bash
set -euo pipefail

# ---------- Config (adjust if you like) ----------
BASE_DIR_DEFAULT="${BASE_DIR_DEFAULT:-/home/omidard/panGEM_pipeline/ESKAPE}"
OUT_DEFAULT="${OUT_DEFAULT:-/home/omidard/genome_qc_results.csv}"
THREADS_DEFAULT="${THREADS_DEFAULT:-64}"

MAMBA_ROOT="${MAMBA_ROOT:-$HOME/micromamba}"
MICROMAMBA_BIN="${MICROMAMBA_BIN:-$(pwd)/bin/micromamba}"   # local bootstrap
ENV_PY="genomeqc_py"    # pandas/biopython runner
ENV_CKM2="ckm2"         # CheckM2 runtime
CKM2_SRC="${HOME}/.local/share/checkm2-src"
CKM2_DB_ROOT="${CKM2_DB_ROOT:-$HOME/databases}"
CKM2_DB_PATH="${CKM2_DB_PATH:-$CKM2_DB_ROOT/CheckM2_database/uniref100.KO.1.dmnd}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/genome_qc.py"

# ---------- Prompts ----------
read -rp "Enter the BASE directory containing genome subfolders [default: ${BASE_DIR_DEFAULT}]: " BASE_DIR_INPUT || true
BASE_DIR="${BASE_DIR_INPUT:-$BASE_DIR_DEFAULT}"
read -rp "Enter the FULL path (CSV/TSV) to save results [default: ${OUT_DEFAULT}]: " OUT_INPUT || true
OUT="${OUT_INPUT:-$OUT_DEFAULT}"
read -rp "Threads to use [default: ${THREADS_DEFAULT}]: " THREADS_INPUT || true
THREADS="${THREADS_INPUT:-$THREADS_DEFAULT}"

# ---------- Helpers ----------
ts() { date '+%m/%d/%Y %I:%M:%S %p'; }
log() { echo "[`ts`] $*"; }

# ---------- Micromamba bootstrap ----------
if ! command -v "$MICROMAMBA_BIN" >/dev/null 2>&1; then
  log "Bootstrapping micromamba (user-local)..."
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba >/dev/null
fi
eval "$("$MICROMAMBA_BIN" shell hook -s bash)"

# ---------- Ensure envs (idempotent) ----------
log "Creating/ensuring env '${ENV_PY}' (pandas/BioPython)..."
"$MICROMAMBA_BIN" create -y -p "${MAMBA_ROOT}/envs/${ENV_PY}" -c conda-forge -c bioconda \
  python=3.11 pandas biopython tqdm >/dev/null || true

log "Cloning CheckM2 (if needed) and creating env '${ENV_CKM2}'..."
if [ ! -d "$CKM2_SRC/.git" ]; then
  git clone --recursive https://github.com/chklovski/checkm2.git "$CKM2_SRC" >/dev/null
fi
# Create env if missing (idempotent)
"$MICROMAMBA_BIN" create -y -p "${MAMBA_ROOT}/envs/${ENV_CKM2}" -c conda-forge -c bioconda \
  "python>3.12" scikit-learn=1.6.1 numpy diamond=2.1.11 tensorflow=2.17 lightgbm pandas scipy \
  prodigal=2.6.3 setuptools requests packaging tqdm >/dev/null || true

# Install/upgrade CheckM2 in that env (idempotent)
"$MICROMAMBA_BIN" run -p "${MAMBA_ROOT}/envs/${ENV_CKM2}" python -m pip install -e "$CKM2_SRC" >/dev/null || true

# ---------- Database ----------
log "Checking CheckM2 database..."
if [ ! -f "$CKM2_DB_PATH" ]; then
  log "DB missing; downloading (~GBs)..."
  "$MICROMAMBA_BIN" run -p "${MAMBA_ROOT}/envs/${ENV_CKM2}" checkm2 database --download --path "$CKM2_DB_ROOT"
else
  echo "      DB OK at $CKM2_DB_PATH"
fi

# ---------- Write the Python (full content, no placeholders) ----------
log "Writing ${PY} ..."
cat > "$PY" <<'PYCODE'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, json, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def parse_args():
    p = argparse.ArgumentParser("Genome QC (GBFF-only) + CheckM2")
    p.add_argument("--base_dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--runner", default="micromamba")
    p.add_argument("--env", default="ckm2")
    p.add_argument("--checkm2_db", default=None)
    p.add_argument("--keep_temp", action="store_true")
    p.add_argument("--reuse_checkm2", action="store_true")
    return p.parse_args()

def find_gbff_in_folder(folder: Path) -> Optional[Path]:
    cands = list(folder.glob("*.gbff")) + list(folder.glob("*.gbk")) + list(folder.glob("*.genbank"))
    if len(cands) == 1:
        return cands[0]
    if len(cands) > 1:
        gbff = [g for g in cands if g.suffix.lower()==".gbff"]
        if len(gbff)==1: return gbff[0]
    return None

def n50(lengths: List[int]) -> int:
    if not lengths: return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths); acc = 0
    for L in lengths:
        acc += L
        if acc >= total/2: return L
    return lengths[-1]

def gc_percent(records) -> float:
    gc = 0; tot = 0
    for r in records:
        s = r.seq.upper()
        tot += len(s); gc += sum(1 for b in s if b in ("G","C"))
    return 0.0 if tot==0 else 100.0*gc/tot

def parse_gbff_and_write_fa(genome: str, gbff: Path, out_fa: Path) -> Tuple[str,int,int,float,Path]:
    recs = list(SeqIO.parse(str(gbff), "genbank"))
    lens = [len(r.seq) for r in recs]
    with out_fa.open("w") as fh:
        for r in recs: SeqIO.write(r, fh, "fasta")
    return genome, len(lens), n50(lens), gc_percent(recs), out_fa

def _norm(x:str)->str:
    stem = Path(x).name
    stem = re.sub(r'\.(fa|fna|fasta|gz)$','',stem,flags=re.I)
    stem = re.sub(r'^\d{1,6}__','',stem)
    return stem

def run_checkm2(fasta_dir:Path,out_dir:Path,threads:int,db:str,runner:str,env:str):
    cmd=[runner,"run","-n",env,"checkm2","predict","--threads",str(threads),
         "--input",str(fasta_dir),"--output-directory",str(out_dir),"--extension","fasta"]
    if db: cmd+=["--database_path",db]
    print("[INFO] Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

def load_ck(out_dir:Path, work_dir:Path)->pd.DataFrame:
    tsv=None
    for c in (out_dir/"quality_report.tsv", out_dir/"quality_report.txt"):
        if c.exists(): tsv=c; break
    if tsv is None: raise FileNotFoundError("No quality_report in "+str(out_dir))
    ck = pd.read_csv(tsv, sep="\t")
    (work_dir/"quality_raw.tsv").write_text(ck.to_csv(sep="\t", index=False))
    cols = {c.lower():c for c in ck.columns}
    name_col = next((cols[c] for c in ["name","bin","bin id","bin_id","genome","filename"] if c in cols), ck.columns[0])
    comp_col = next((c for c in ck.columns if "completeness" in c.lower() or c.lower()=="complete"), None)
    cont_col = next((c for c in ck.columns if "contamination" in c.lower() or "contam"==c.lower()), None)
    if comp_col is None or cont_col is None:
        raise RuntimeError("Cannot find completeness/contamination columns in CheckM2 output")
    nm_path = work_dir/"name_map.json"
    name_map:Dict[str,str] = json.loads(nm_path.read_text())
    rev = {_norm(pref):genome for pref,genome in name_map.items()}
    ck["_key"] = ck[name_col].astype(str).map(_norm)
    ck["genome"] = ck["_key"].map(rev)
    unmapped = ck[ck["genome"].isna()][name_col].dropna().tolist()
    if unmapped:
        (work_dir/"unmapped_bins.txt").write_text("\n".join(map(str,unmapped))+"\n")
        print(f"[WARN] {len(unmapped)} CheckM2 bins unmapped -> {work_dir/'unmapped_bins.txt'}", flush=True)
    out = ck.loc[ck["genome"].notna(), ["genome", comp_col, cont_col]].copy()
    out.rename(columns={comp_col:"completeness", cont_col:"contamination"}, inplace=True)
    out["completeness"]=pd.to_numeric(out["completeness"],errors="coerce")
    out["contamination"]=pd.to_numeric(out["contamination"],errors="coerce")
    return out

def main():
    a = parse_args()
    base = Path(a.base_dir).resolve()
    out = Path(a.out).resolve()
    thr = max(1,int(a.threads))
    runner = str(Path(a.runner).resolve() if "/" in a.runner else a.runner)
    work = Path(tempfile.mkdtemp(prefix="genome_qc_"))
    fasta_dir = work/"assemblies"; fasta_dir.mkdir(parents=True, exist_ok=True)
    ck_out = work/"checkm2"; ck_out.mkdir(parents=True, exist_ok=True)

    print("[INFO] Scanning for genome folders with GBFF ...", flush=True)
    genomes=[]; gbff_by={}
    for sub in sorted(base.iterdir()):
        if not sub.is_dir(): continue
        gbff = find_gbff_in_folder(sub)
        if gbff is None:
            print(f"[WARN] No GBFF in: {sub}", flush=True)
            continue
        genomes.append(sub.name); gbff_by[sub.name]=gbff
    if not genomes:
        print("[FATAL] No genomes with GBFF found", flush=True); sys.exit(2)
    print(f"[INFO] Found {len(genomes)} genome folders with GBFF.", flush=True)

    # Export FASTAs + per-genome stats (multiproc) with stable prefixes
    name_map={}
    jobs=[]
    with ProcessPoolExecutor(max_workers=thr) as ex:
        for i,g in enumerate(genomes, start=1):
            pref=f"{i:05d}__{g}"; name_map[pref]=g
            jobs.append(ex.submit(parse_gbff_and_write_fa, g, gbff_by[g], fasta_dir/f"{pref}.fasta"))
        rows=[]
        for fut in tqdm(as_completed(jobs), total=len(jobs), desc="Genomes"):
            g,nc,n50,gc,fa = fut.result()
            rows.append({"genome":g,"num_contigs":int(nc),"N50":int(n50),
                         "gc_percent":round(float(gc),3),
                         "gbff_path":str(gbff_by[g]),"fasta_path":str(fa)})
    (work/"name_map.json").write_text(json.dumps(name_map))
    df = pd.DataFrame(rows).sort_values("genome").reset_index(drop=True)

    # Run CheckM2
    try:
        run_checkm2(fasta_dir, ck_out, thr, a.checkm2_db, runner, a.env)
        ck = load_ck(ck_out, work)
        df = df.merge(ck, on="genome", how="left")
    except Exception as e:
        print(f"[WARN] CheckM2 step failed or not merged: {e}", flush=True)
        if "completeness" not in df: df["completeness"]=pd.NA
        if "contamination" not in df: df["contamination"]=pd.NA

    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower()==".tsv": df.to_csv(out, sep="\t", index=False)
    else: df.to_csv(out, index=False)

    print(f"[OK] Saved: {out}", flush=True)
    print(f"[STATS] genomes_detected={len(genomes)} rows_written={len(df)}", flush=True)

    if not a.keep_temp:
        print(f"[INFO] Cleaning temp at {work}", flush=True)
        shutil.rmtree(work, ignore_errors=True)
    else:
        print(f"[INFO] Temp kept at: {work}", flush=True)

if __name__=="__main__":
    main()
PYCODE
chmod +x "$PY"
# sanity: ensure file is not a tiny placeholder
MIN_SIZE=5000
ACTUAL_SIZE=$(wc -c < "$PY" || echo 0)
if [ "$ACTUAL_SIZE" -lt "$MIN_SIZE" ]; then
  echo "[FATAL] genome_qc.py looks too small (${ACTUAL_SIZE} bytes). Aborting."
  exit 3
fi

# ---------- Run ----------
log "Running genome_qc.py using modern env..."
set +e
OUT_LOG="$("$MICROMAMBA_BIN" run -p "${MAMBA_ROOT}/envs/${ENV_PY}" python "$PY" \
  --base_dir "$BASE_DIR" \
  --out "$OUT" \
  --threads "$THREADS" \
  --runner "$MICROMAMBA_BIN" \
  --env "${ENV_CKM2}" \
  --checkm2_db "$CKM2_DB_PATH" \
  2>&1)"
RC=$?
set -e
echo "$OUT_LOG"
if [ $RC -ne 0 ]; then
  echo "[FATAL] genome_qc.py failed with exit code $RC"
  exit $RC
fi

# ---------- Sanity check the output ----------
if [ ! -s "$OUT" ]; then
  echo "[FATAL] Output file not created: $OUT"
  exit 4
fi

# Expect at least 10 rows (1 header + 9 data) — adjust if needed
ROWS=$(wc -l < "$OUT" | awk '{print $1}')
if [ "$ROWS" -lt 10 ]; then
  echo "[FATAL] Output looks incomplete (only $((ROWS-1)) data rows)."
  echo "        See the run log above for [STATS] and any [WARN]/[FATAL] lines."
  exit 5
fi

echo
echo "[DONE] QC table written to: $OUT"
