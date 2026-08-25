#!/usr/bin/env bash
#
# Smoke test for the minimal example dataset.
#
# Builds a compact lnc-seeker cache from the three bundled GSE235167
# samples across the six Supplementary Note 6 loci, then verifies the
# outputs. This confirms the tool is installed, configured and able to
# read the bundled BAM + GTF data end to end.
#
# Usage:  ./run_example.sh
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# --- Locate the lnc-seeker-collect binary ---
BIN="${LNC_SEEKER_COLLECT:-lnc-seeker-collect}"
if ! command -v "$BIN" >/dev/null 2>&1; then
    REPO_ROOT="$(cd "$HERE/../../" && pwd)"
    for cand in "$REPO_ROOT/target/debug/lnc-seeker-collect" \
                "$REPO_ROOT/target/release/lnc-seeker-collect"; do
        if [ -x "$cand" ]; then BIN="$cand"; break; fi
    done
fi
if ! command -v "$BIN" >/dev/null 2>&1 && [ ! -x "${BIN:-}" ]; then
    echo "ERROR: lnc-seeker-collect not found." >&2
    echo "Build it first with:  cargo build -p lnc-seeker-collect" >&2
    echo "(or set LNC_SEEKER_COLLECT=/path/to/lnc-seeker-collect)" >&2
    exit 1
fi
echo "Using binary: $BIN"

# --- Run collection ---
rm -rf collect_output
"$BIN" --config config_example.cfg

# --- Verify outputs ---
echo
echo "=== collect_output/ ==="
ls -lh collect_output

fail=0
if [ -s collect_output/dictionary.json ]; then
    echo "OK: dictionary.json present and non-empty"
else
    echo "FAIL: dictionary.json missing or empty"; fail=1
fi

# Collection writes one compact cache per locus: <cohort>_<gene>_<hash>.lnc_cache.bin
n_cache="$(find collect_output -name '*.lnc_cache.bin' 2>/dev/null | wc -l)"
if [ "$n_cache" -ge 6 ]; then
    echo "OK: $n_cache per-locus cache files present (expected one per locus)"
else
    echo "FAIL: expected >=6 .lnc_cache.bin files, found $n_cache"; fail=1
fi

# Built-in sanity check: PCA3 is prostate-specific, so in these breast
# samples its cache should be tiny (a handful of reads) relative to GAPDH.
echo
echo "PCA3 vs GAPDH cache size (tissue-specificity sanity check):"
ls -lh collect_output/*_PCA3_*.lnc_cache.bin collect_output/*_GAPDH_*.lnc_cache.bin 2>/dev/null \
    | awk '{print "    " $NF "  (" $5 ")"}'

echo
if [ "$fail" -eq 0 ]; then
    echo "SUCCESS: minimal dataset ran end to end."
    echo "Next: point a lnc-seeker config.json at collect_output/ to open it in the app."
else
    echo "FAILURE: see messages above." >&2
fi
exit "$fail"
