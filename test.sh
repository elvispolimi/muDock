#!/usr/bin/env bash
#
# Runs every ligand/protein pair found in the muDock test sets through both muDock's X-Score
# (x_score_bench) and the original XScore, writing a side-by-side comparison of the raw terms that are
# currently implemented (VDW HB HP RT) plus the affinity they predict (PKD) to output.txt.
#
# PKD is XScore's pkd1, the HPScore prediction -- of the three -log(Kd) values XScore reports, the only
# one whose regression uses exactly the four terms muDock implements:
#     pKd = 3.441 + 0.004*VDW + 0.054*HB + 0.009*HP - 0.061*RT
# XScore's other two (HMScore, HSScore) need the HM / HS terms, which muDock does not compute, so they
# are deliberately not compared.
#
# On top of the human-readable comparison, every pair is checked term by term: both values are ROUNDED
# to three decimals and the rounded results must be equal. Rounding, not truncation -- 1.9325 becomes
# 1.933, not 1.932. Note the consequence: two values may agree to within 1e-4 and still round to
# different results when they straddle a x.xxx5 boundary, so a MISMATCH is not by itself evidence of a
# larger disagreement. The raw delta is printed next to every mismatch so that case is recognisable.
#
# Terms are being re-introduced one at a time; extend the *_RE regexes below as new terms come online.
#
# Usage:  ./test.sh          ; exit status is 0 only if every pair matches
#
set -u

# printf must use '.' as the decimal separator for the rounding below to behave.
export LC_ALL=C

# Number of decimals the comparison rounds to.
ROUND_DECIMALS=3

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- locations (edit here if your layout differs) ---
PROT_DIR="$HOME/xscore-dev/srcFiles/proteins/mudock_proteins"
LIG_DIR="$HOME/xscore-dev/srcFiles/ligands/mudock_ligands"
XTOOL_DIR="$HOME/xscore-dev/srcFiles/xtool"           # holds the score_<name>.input files
BENCH="$ROOT/build/application/bench/x_score_bench"
XSCORE="$HOME/xscore-dev/xscore_v1.3/c++/xscore"
export XSCORE_PARAMETER="$HOME/xscore-dev/xscore_v1.3/parameter/"

OUT="$ROOT/output.txt"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Extract the currently-implemented terms (VDW, HB, HP, RT) and the resulting affinity (PKD) from a
# tool's output. Each value is printed by muDock and XScore as "VDW=<value>" / "HB=<value>" /
# "HP=<value>" / "RT=<value>" / "PKD=<value>"; we pull them individually so the two tools can print
# different sets of terms (XScore also emits HM/HS) without breaking the comparison.
terms_of() {
  # $1 = raw multi-line tool output ; echoes "VDW=<> HB=<> HP=<> RT=<> PKD=<>"
  local text="$1" vdw hb hp rt pkd
  vdw="$(printf '%s\n' "$text" | grep -oE 'VDW=[^ ]+' | head -1)"
  hb="$(printf '%s\n' "$text"  | grep -oE 'HB=[^ ]+'  | head -1)"
  hp="$(printf '%s\n' "$text"  | grep -oE 'HP=[^ ]+'  | head -1)"
  rt="$(printf '%s\n' "$text"  | grep -oE 'RT=[^ ]+'  | head -1)"
  pkd="$(printf '%s\n' "$text" | grep -oE 'PKD=[^ ]+' | head -1)"
  echo "${vdw:-VDW=?} ${hb:-HB=?} ${hp:-HP=?} ${rt:-RT=?} ${pkd:-PKD=?}"
}

# Pull a single term's value out of a "VDW=<> HB=<> HP=<> RT=<> PKD=<>" string.
value_of() {
  # $1 = terms string, $2 = term name ; echoes the bare value ("?" when the tool printed nothing)
  printf '%s\n' "$1" | tr ' ' '\n' | sed -n "s/^$2=//p" | head -1
}

# Round a value to ROUND_DECIMALS decimals. A missing value stays "?" so it can never compare equal.
round() {
  # $1 = value
  case "$1" in
    ''|'?') printf '?' ;;
    *) printf "%.${ROUND_DECIMALS}f" "$1" ;;
  esac
}

# --- sanity checks ---
if [ ! -x "$BENCH" ]; then
  echo "x_score_bench not found at $BENCH — building it ..."
  cmake --build "$ROOT/build" --target x_score_bench -j4 || { echo "build failed"; exit 1; }
fi
[ -x "$XSCORE" ] || { echo "xscore binary not found at $XSCORE"; exit 1; }
[ -d "$PROT_DIR" ] || { echo "protein dir not found: $PROT_DIR"; exit 1; }
[ -d "$LIG_DIR" ]  || { echo "ligand dir not found: $LIG_DIR"; exit 1; }

: > "$OUT"

run_xscore() {
  # $1 = name, $2 = protein pdb, $3 = ligand mol2 ; echoes the full xscore output
  local name="$1" pdb="$2" lig="$3" inp
  if [ -f "$XTOOL_DIR/score_${name}.input" ]; then
    # use the validated input (relative paths -> run from its directory)
    ( cd "$XTOOL_DIR" && "$XSCORE" "score_${name}.input" ) 2>&1
  else
    # generate a minimal input for this pair (defaults for the coefficients)
    inp="$TMP/${name}.input"
    cat > "$inp" <<EOF
FUNCTION            SCORE
RECEPTOR_PDB_FILE   $pdb
LIGAND_MOL2_FILE    $lig
OUTPUT_TABLE_FILE   $TMP/${name}.table
OUTPUT_LOG_FILE     $TMP/${name}.log
NUMBER_OF_HITS      1
HITS_DIRECTORY      $TMP/${name}.mdb
SHOW_ATOM_BIND_SCORE    NO
APPLY_HPSCORE       YES
APPLY_HMSCORE       YES
APPLY_HSSCORE       YES
APPLY_CHEMICAL_RULES    NO
EOF
    "$XSCORE" "$inp" 2>&1
  fi
}

count=0
ok_count=0
fail_count=0
failed_pairs=""
for pdb in "$PROT_DIR"/*_protein.pdb; do
  [ -e "$pdb" ] || continue
  name="$(basename "$pdb" _protein.pdb)"
  lig="$LIG_DIR/${name}_ligand.mol2"
  if [ ! -f "$lig" ]; then
    echo "skip $name (no ligand at $lig)"
    continue
  fi

  echo "scoring $name ..."
  md="$(terms_of "$("$BENCH" -p "$pdb" -l "$lig" 2>&1)")"
  xs="$(terms_of "$(run_xscore "$name" "$pdb" "$lig")")"

  # Term-by-term check: round both values to ROUND_DECIMALS decimals and require equality.
  pair_ok=1
  mismatches=""
  bad_terms=""
  for term in VDW HB HP RT PKD; do
    mv="$(value_of "$md" "$term")"
    xv="$(value_of "$xs" "$term")"
    mr="$(round "$mv")"
    xr="$(round "$xv")"
    if [ "$mr" != '?' ] && [ "$mr" = "$xr" ]; then
      continue
    fi
    pair_ok=0
    bad_terms="$bad_terms$term,"
    if [ "$mr" = '?' ] || [ "$xr" = '?' ]; then
      mismatches="$mismatches
    MISMATCH  $term  mudock=${mv:-?} xscore=${xv:-?}  (term missing from one of the two tools)"
    else
      delta="$(awk -v a="$mv" -v b="$xv" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.3e", d}')"
      mismatches="$mismatches
    MISMATCH  $term  mudock=$mv -> $mr   xscore=$xv -> $xr   (raw delta $delta)"
    fi
  done

  {
    echo "$name"
    echo "    mudock    ${md:-<no output>}"
    echo "    xscore    ${xs:-<no output>}"
    if [ "$pair_ok" -eq 1 ]; then
      echo "    check     OK (all terms equal when rounded to $ROUND_DECIMALS decimals)"
    else
      printf '%s\n' "    check     FAILED${mismatches}"
    fi
    echo
  } >> "$OUT"

  if [ "$pair_ok" -eq 1 ]; then
    ok_count=$((ok_count + 1))
  else
    fail_count=$((fail_count + 1))
    failed_pairs="$failed_pairs $name(${bad_terms%,})"
  fi
  count=$((count + 1))
done

# --- summary, written to both the report and the terminal ---
summary="$(
  echo "=== SUMMARY ==="
  echo "pairs compared : $count"
  echo "matching       : $ok_count"
  echo "mismatching    : $fail_count"
  echo "criterion      : every term equal after rounding to $ROUND_DECIMALS decimals"
  if [ "$fail_count" -gt 0 ]; then
    echo "mismatching pairs:$failed_pairs"
  fi
)"
printf '%s\n' "$summary" >> "$OUT"
printf '%s\n' "$summary"

echo "Done: $count pair(s) written to $OUT"

if [ "$fail_count" -gt 0 ]; then
  exit 1
fi
exit 0