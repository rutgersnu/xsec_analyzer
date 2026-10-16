#!/usr/bin/env bash
# nuNCQE chain: annie_ncqe_prep (parallel, one job per MC file) -> hadd -> univmake -> ncqe_xsec
#
# Usage: ncqe/run_ncqe.sh TAG MC_DIR MC_POT SYSTCALC_CONF [options]
#   TAG            name of this run (output subdirectory)
#   MC_DIR         directory with ToolAnalysis MC files MC_<n>.root
#   MC_POT         simulated POT of the files that are processed (e.g. 4.41733e20 for all 3922 MC_GENIE-O16 files)
#   SYSTCALC_CONF  e.g. ncqe/systcalc_cv.conf (no unisims) or ncqe/systcalc_44knob.conf
# Options:
#   --unisims        also write the 8 unisim slots (44-knob reprocessing)
#   --nfiles N       only the first N files (by file number); MC_POT must then be the POT of that subset
#   --list FILE      explicit list of file names in MC_DIR (one per line) instead of all MC_*.root
#   --jobs J         parallel prep jobs (default: nproc)
#   --inputs FILE    data-side inputs (default: ncqe/ncqe_inputs.conf)
#   --out DIR        output base directory (default: $NCQE_OUT or ~/builds/ncqe_xsa); results go to DIR/TAG
#   --data-pot POT   beam-on POT (default 2.572970591694e20, FY22-23)
#   --data-trigs N   beam-on triggers (default 60624559)
#   --skip-prep      reuse the merged MC file of a previous run with the same TAG
set -euo pipefail

if [ $# -lt 4 ]; then sed -n '2,20p' "$0"; exit 1; fi
REPO=$(cd "$(dirname "$0")/.." && pwd)
TAG=$1; MC_DIR=$(cd "$2" && pwd); MC_POT=$3; SYSTCALC=$(realpath "$4"); shift 4

UNISIM_FLAG=""; NFILES=""; LIST=""; JOBS=$(nproc); INPUTS="$REPO/ncqe/ncqe_inputs.conf"
OUTBASE=${NCQE_OUT:-$HOME/builds/ncqe_xsa}; DATA_POT=2.572970591694e20; DATA_TRIGS=60624559; SKIP_PREP=0
while [ $# -gt 0 ]; do
  case $1 in
    --unisims) UNISIM_FLAG="--unisims"; shift ;;
    --nfiles) NFILES=$2; shift 2 ;;
    --list) LIST=$(realpath "$2"); shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;;
    --inputs) INPUTS=$(realpath "$2"); shift 2 ;;
    --out) OUTBASE=$2; shift 2 ;;
    --data-pot) DATA_POT=$2; shift 2 ;;
    --data-trigs) DATA_TRIGS=$2; shift 2 ;;
    --skip-prep) SKIP_PREP=1; shift ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done

export STV_ANALYSIS_DIR=$REPO
OUT=$OUTBASE/$TAG
mkdir -p "$OUT/prep"
MERGED=$OUT/ncqe_mc_$TAG.root
make -C "$REPO" ncqe > "$OUT/build.log" 2>&1 || { cat "$OUT/build.log"; exit 1; }

# ------------------------------------------------------------------ 1. file list and per-file prep
if [ "$SKIP_PREP" = 0 ]; then
  if [ -n "$LIST" ]; then
    grep -v '^\s*$' "$LIST" | grep -v '^#' > "$OUT/mc_files.txt"
  else
    ls "$MC_DIR" | grep -E '^MC_[0-9]+\.root$' | sort -t_ -k2 -n > "$OUT/mc_files.txt"
  fi
  if [ -n "$NFILES" ]; then head -n "$NFILES" "$OUT/mc_files.txt" > "$OUT/mc_files.tmp"; mv "$OUT/mc_files.tmp" "$OUT/mc_files.txt"; fi
  N_LISTED=$(wc -l < "$OUT/mc_files.txt")
  echo "[$TAG] prep of $N_LISTED files from $MC_DIR with $JOBS jobs"
  rm -f "$OUT"/prep/*.ncqe.root

  prep_one() {
    local name=$1
    "$REPO/bin/annie_ncqe_prep" "$MC_DIR/$name" "$OUT/prep/${name%.root}.ncqe.root" $UNISIM_FLAG
  }
  export -f prep_one; export REPO MC_DIR OUT UNISIM_FLAG
  xargs -a "$OUT/mc_files.txt" -P "$JOBS" -n 1 bash -c 'prep_one "$1"' _ > "$OUT/prep.log" 2>&1 || true

  N_DONE=$(grep -c '^PREP ' "$OUT/prep.log" || true)
  N_SKIP=$(grep -c '^SKIP ' "$OUT/prep.log" || true)
  echo "[$TAG] prepped $N_DONE, skipped $N_SKIP (unreadable), of $N_LISTED"
  grep -E '^(SKIP|WARNING|terminate|.*rror)' "$OUT/prep.log" | head -20 || true
  if [ $((N_DONE + N_SKIP)) -ne "$N_LISTED" ]; then echo "[$TAG] prep failures, see $OUT/prep.log"; exit 1; fi
  # Totals over files (columns of the PREP summary lines)
  awk '/^PREP /{for(i=3;i<NF;i+=2) s[$i]+=$(i+1)} END{printf "[totals]"; for(k in s) printf " %s %d", k, s[k]; printf "\n"}' \
    "$OUT/prep.log" | tee "$OUT/prep_totals.txt"

  find "$OUT/prep" -name '*.ncqe.root' | sort > "$OUT/prep_outputs.txt"
  # shellcheck disable=SC2046
  hadd -f -k "$MERGED" $(cat "$OUT/prep_outputs.txt") > "$OUT/hadd.log" 2>&1 || { tail "$OUT/hadd.log"; exit 1; }
fi

# ------------------------------------------------------------------ 2. placeholders and framework bookkeeping
"$REPO/bin/annie_ncqe_prep" --placeholder "$OUT/ncqe_onbnb_placeholder.root" > /dev/null
"$REPO/bin/annie_ncqe_prep" --placeholder "$OUT/ncqe_extbnb_placeholder.root" > /dev/null
FILES=$OUT/files_$TAG.txt
FP=$OUT/file_properties_$TAG.txt
printf '%s\n' "$OUT/ncqe_onbnb_placeholder.root" "$OUT/ncqe_extbnb_placeholder.root" "$MERGED" > "$FILES"
{
  echo "# nuNCQE run $TAG: placeholders carry the beam-on POT and triggers; N_obs, D_CIT, D_sky come from $INPUTS"
  echo "$OUT/ncqe_onbnb_placeholder.root 1 onBNB $DATA_TRIGS $DATA_POT"
  echo "$OUT/ncqe_extbnb_placeholder.root 1 extBNB 1 0"
  echo "$MERGED 1 numuMC $MC_POT"
} > "$FP"

# ------------------------------------------------------------------ 3. universes (fresh file: cached totals are reused silently)
UNIV=$OUT/univmake_$TAG.root
rm -f "$UNIV"
echo "[$TAG] univmake"
"$REPO/bin/univmake" "$FILES" "$REPO/ncqe/bin_config_ncqe.txt" "$UNIV" "$FP" > "$OUT/univmake.log" 2>&1 \
  || { tail -20 "$OUT/univmake.log"; exit 1; }

# ------------------------------------------------------------------ 4. extraction and uncertainties
echo "[$TAG] ncqe_xsec"
"$REPO/bin/ncqe_xsec" "$UNIV" "$SYSTCALC" "$INPUTS" "$FP" "$OUT/ncqe_xsec_$TAG.json" > "$OUT/ncqe_xsec_full.log" 2>&1 \
  || { tail -30 "$OUT/ncqe_xsec_full.log"; exit 1; }
# The framework prints many DEBUG lines; keep the summary separately
tr '\r' '\n' < "$OUT/ncqe_xsec_full.log" | grep -v -E 'DEBUG|^\s*$|PROCESSING' > "$OUT/ncqe_xsec_$TAG.log"
sed -n '/==== nuNCQE extraction/,$p' "$OUT/ncqe_xsec_$TAG.log"
echo "[$TAG] outputs in $OUT"
