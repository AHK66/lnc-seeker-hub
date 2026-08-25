#!/usr/bin/env python3
r"""Reproducibility of tumor-specific splice junctions across independent cohorts.

Answers the reviewer question on the lncRNA Seeker Hub splice-junction
evaluation (Supplementary Note 6): how many tumor-specific (and novel) splice
junctions found in the TCGA-BRCA tumor cohort are reproducibly detected in an
independent breast-cancer cohort (GSE235167 PDX)?

Inputs are the per-cohort, per-locus J_LSH junction TSVs already produced by the
Note-6 pipeline (run_comparison.sh), so no BAM re-streaming is required. The
junction definition is reused verbatim from compare_junctions.py (normalize
chrom -> region-filter fully-inside -> threshold on cohort-summed reads); the
script self-checks against the published Note-6 totals before reporting.

Definitions (junction identity = unstranded (chrom, start, end)):
    tumor-specific      = J_t(BRCA tumor) - J_t(BRCA normal)
    breast-tumor-spec.  = J_t(BRCA tumor) - (J_t(BRCA normal) | J_t(PRAD))
                          (absent in BOTH matched normal AND prostate ->
                           excludes broadly/housekeeping-expressed junctions)
    reproduced          = above  |  J_t(PDX)          [independent breast cohort]
    prostate overlap    = tumor-specific  &  J_t(PRAD) [unrelated-tissue control]
    novel               = above  -  (RefSeq | GENCODE lncRNA)

Usage:
    python3 tumor_specific_reproducibility.py [--base PATH] [--out PATH]

  --base  repo root holding the Note-6 pipeline outputs
         (junction_comparison_output{,_brca_normal,_prad,_gse235167}/)
  --out   directory for the result CSVs (default: tumor_specific_repro/ here)
"""

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compare_junctions import load_junctions, filter_by_region, build_sets  # noqa: E402

# Default: the local checkout where the Note-6 pipeline was run and its
# per-region J_LSH TSVs + resolved config + annotation TSVs live.
DEFAULT_BASE = "/path/to/lnc-seeker-hub"
GENES = ["ACTB", "PCA3", "CD44", "MALAT1", "GAPDH", "TP53"]
THRESHOLDS = [1, 2, 5, 10]
# Note-6 published |J_LSH| totals per threshold (Supplementary Table 6.3),
# used only to verify the junction definition matches the pipeline exactly.
NOTE6 = {"tumor": {1: 323, 2: 272, 5: 216, 10: 181},
         "normal": {1: 245, 2: 191, 5: 152, 10: 128},
         "prad": {1: 320, 2: 259, 5: 198, 10: 174},
         "pdx": {1: 806, 2: 533, 5: 255, 10: 185}}


def cohort_dirs(base: str) -> dict:
    return {
        "tumor":  f"{base}/junction_comparison_output",
        "normal": f"{base}/junction_comparison_output_brca_normal",
        "pdx":    f"{base}/junction_comparison_output_gse235167",
        "prad":   f"{base}/junction_comparison_output_prad",
    }


def load_regions(base: str) -> dict:
    cfg = json.load(open(f"{base}/junction_comparison_output/config_resolved.json"))
    return {r["gene"]: f'{r["chrom"]}:{r["start"]}-{r["end"]}' for r in cfg["resolved_regions"]}


def lnc_set(directories: dict, regions: dict, cohort: str, gene: str, threshold: int) -> set:
    path = f'{directories[cohort]}/lncseeker/per_region/{gene}_junctions.tsv'
    js = filter_by_region(load_junctions(path), regions[gene])
    return build_sets(js, threshold, stranded=False)


def load_known(directories: dict) -> set:
    known = set()
    for tag in ("refseq_junctions.tsv", "gencode_junctions.tsv"):
        p = f'{directories["tumor"]}/stringtie/{tag}'
        for c, s, e, _strand, _n in load_junctions(p):
            known.add((c, s, e))
    return known


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE, help="repo root with Note-6 pipeline outputs")
    ap.add_argument("--out", default=os.path.join(HERE, "tumor_specific_repro"), help="output CSV directory")
    args = ap.parse_args()

    dirs = cohort_dirs(args.base)
    regions = load_regions(args.base)
    known = load_known(dirs)
    os.makedirs(args.out, exist_ok=True)

    per_locus = []
    baseline = {c: {t: 0 for t in THRESHOLDS} for c in dirs}

    for gene in GENES:
        for t in THRESHOLDS:
            tumor  = lnc_set(dirs, regions, "tumor", gene, t)
            normal = lnc_set(dirs, regions, "normal", gene, t)
            pdx    = lnc_set(dirs, regions, "pdx", gene, t)
            prad   = lnc_set(dirs, regions, "prad", gene, t)
            for cohort, s in (("tumor", tumor), ("normal", normal), ("pdx", pdx), ("prad", prad)):
                baseline[cohort][t] += len(s)

            ts = tumor - normal                 # tumor-specific (loose)
            ts_strict = ts - prad               # breast-tumor-specific
            neg = ts & prad                     # prostate overlap
            repro = ts & pdx
            repro_strict = ts_strict & pdx
            novel_ts = {j for j in ts if j not in known}
            novel_repro = novel_ts & pdx

            per_locus.append({
                "gene": gene, "threshold": t,
                "tumor": len(tumor), "normal": len(normal),
                "tumor_specific": len(ts),
                "neg_in_prad": len(neg),
                "breast_tumor_specific": len(ts_strict),
                "repro_in_pdx": len(repro),
                "repro_frac": round(len(repro) / len(ts), 4) if ts else 0.0,
                "breast_ts_repro_in_pdx": len(repro_strict),
                "breast_ts_repro_frac": round(len(repro_strict) / len(ts_strict), 4) if ts_strict else 0.0,
                "pdx": len(pdx), "prad": len(prad),
                "novel_ts": len(novel_ts),
                "novel_repro_in_pdx": len(novel_repro),
                "novel_repro_frac": round(len(novel_repro) / len(novel_ts), 4) if novel_ts else 0.0,
            })

    with open(os.path.join(args.out, "tumor_specific_per_locus.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_locus[0].keys()))
        w.writeheader()
        w.writerows(per_locus)

    overall = []
    for t in THRESHOLDS:
        rows = [r for r in per_locus if r["threshold"] == t]
        S = lambda k: sum(r[k] for r in rows)
        ts, bts = S("tumor_specific"), S("breast_tumor_specific")
        rep, brep = S("repro_in_pdx"), S("breast_ts_repro_in_pdx")
        nt, nr = S("novel_ts"), S("novel_repro_in_pdx")
        overall.append({
            "threshold": t,
            "tumor_specific": ts,
            "repro_in_pdx": rep,
            "repro_frac": round(rep / ts, 4) if ts else 0.0,
            "neg_in_prad": S("neg_in_prad"),
            "breast_tumor_specific": bts,
            "breast_ts_repro_in_pdx": brep,
            "breast_ts_repro_frac": round(brep / bts, 4) if bts else 0.0,
            "novel_ts": nt,
            "novel_repro_in_pdx": nr,
            "novel_repro_frac": round(nr / nt, 4) if nt else 0.0,
        })
    with open(os.path.join(args.out, "tumor_specific_overall.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(overall[0].keys()))
        w.writeheader()
        w.writerows(overall)

    # ---- self-check against Note-6 published totals ----
    ok_all = True
    print("=== BASELINE CHECK: |J_LSH| over 6 loci (vs Note-6 Table 6.3) ===")
    for c in ("tumor", "normal", "prad", "pdx"):
        vals = "  ".join(f"{baseline[c][t]:<4}" for t in THRESHOLDS)
        match = all(baseline[c][t] == NOTE6[c][t] for t in THRESHOLDS)
        ok_all &= match
        print(f"{c:8} {vals}   {'OK' if match else 'MISMATCH'}")
    if not ok_all:
        print("WARNING: baseline does NOT match Note-6; do not trust the results below.", file=sys.stderr)
        sys.exit(1)
    print("baseline MATCHES Note-6\n")

    print("=== REPRODUCIBILITY: BRCA tumor-specific vs independent PDX breast cohort (unstranded) ===")
    print(f"{'t':>3} | {'tum-spec':>8} {'@PDX':>5} {'%':>6} | {'@PRAD':>5} | "
          f"{'br-tspec':>8} {'@PDX':>5} {'%':>6} | {'novel':>5} {'@PDX':>5} {'%':>6}")
    for r in overall:
        print(f"{r['threshold']:>3} | {r['tumor_specific']:>8} {r['repro_in_pdx']:>5} "
              f"{r['repro_frac']*100:>5.1f}% | {r['neg_in_prad']:>5} | "
              f"{r['breast_tumor_specific']:>8} {r['breast_ts_repro_in_pdx']:>5} "
              f"{r['breast_ts_repro_frac']*100:>5.1f}% | "
              f"{r['novel_ts']:>5} {r['novel_repro_in_pdx']:>5} {r['novel_repro_frac']*100:>5.1f}%")

    print(f"\nWrote: {os.path.join(args.out, 'tumor_specific_per_locus.csv')}")
    print(f"       {os.path.join(args.out, 'tumor_specific_overall.csv')}")


if __name__ == "__main__":
    main()
