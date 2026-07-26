#!/usr/bin/env bash
set -euo pipefail

# Junction comparison pipeline: lnc-seeker vs StringTie (cohort-oriented)
# Usage: ./junction_comparison/run_comparison.sh [config_path] [output_dir]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG="${1:-$REPO_DIR/junction_comparison/config_comparison.cfg}"
OUTDIR="${2:-$(python3 -c "import sys; sys.path.insert(0,'$SCRIPT_DIR'); from parse_config import parse_config; print(parse_config('$CONFIG')['output_dir'])")}"

mkdir -p "$OUTDIR"/{stringtie,lncseeker,comparison,logs}

echo "=== Parsing config and resolving gene regions ==="
# Use gtf_stringtie (uncompressed) for gene region resolution
python3 -c "
import json, sys
sys.path.insert(0, '$SCRIPT_DIR')
from parse_config import parse_config, resolve_gene_regions

config = parse_config('$CONFIG')
gtf_for_resolve = config.get('gtf_stringtie') or config.get('gtf')
config['resolved_regions'] = resolve_gene_regions(gtf_for_resolve, config['gene_regions'])
with open('$OUTDIR/config_resolved.json', 'w') as f:
    json.dump(config, f, indent=2)
" 2>&1 | tee "$OUTDIR/logs/config.log"
python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(json.dumps(c, indent=2))"

# Extract key values from resolved config
STR=$(python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(c.get('stringtie', 'stringtie'))")
GTF=$(python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(c['gtf'])")
GTF_ST=$(python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(c.get('gtf_stringtie', c['gtf']))")
GENOME=$(python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(c['genome'])")
BAMS=$(python3 -c "
import json; c=json.load(open('$OUTDIR/config_resolved.json'))
bams = list(c['bam_files'])
print(' '.join(bams))
")
REGIONS=$(python3 -c "
import json; c=json.load(open('$OUTDIR/config_resolved.json'))
regions = c.get('resolved_regions', [])
for r in regions:
    print(f\"{r['chrom']}:{r['start']}-{r['end']}\")
")

echo "StringTie: $STR"
echo "GTF: $GTF"
echo "Genome: $GENOME"
echo "BAM files: $(echo "$BAMS" | wc -w)"
echo "Regions: $(echo "$REGIONS" | wc -l)"
echo ""

# Step 1: Run StringTie on each BAM independently
echo "=== Step 1: Running StringTie on each BAM ==="
ST_GTF_LIST=""
for bam in $BAMS; do
    sample_name=$(basename "$bam" | sed 's/_Aligned\.sortedByCoord\.out\.bam//;s/\.bam$//')
    st_gtf="$OUTDIR/stringtie/${sample_name}.gtf"
    if [ ! -f "$st_gtf" ]; then
        echo "  StringTie: $sample_name"
        $STR -o "$st_gtf" "$bam" 2>"$OUTDIR/logs/stringtie_${sample_name}.log"
    fi
    ST_GTF_LIST="$ST_GTF_LIST $st_gtf"
done

# Step 2: Merge all StringTie GTFs
echo "=== Step 2: Merging StringTie results ==="
ST_MERGED="$OUTDIR/stringtie/merged.gtf"
if [ ! -f "$ST_MERGED" ]; then
    $STR --merge -o "$ST_MERGED" $ST_GTF_LIST \
        2>"$OUTDIR/logs/stringtie_merge.log"
fi

# Step 3: Extract junctions from merged StringTie GTF
echo "=== Step 3: Extracting StringTie junctions ==="
ST_JUNCTIONS="$OUTDIR/stringtie/junctions.tsv"
python3 "$SCRIPT_DIR/extract_stringtie_junctions.py" "$ST_MERGED" -o "$ST_JUNCTIONS"

# Step 4: (Region filtering is handled by compare_junctions.py --region flag)
echo "=== Step 4: Junctions will be filtered during comparison ==="
ST_JUNCTIONS_FILTERED="$ST_JUNCTIONS"

# Step 5: Run lnc-seeker per-region analysis (indexed BAM, fast)
echo "=== Step 5: Running lnc-seeker per-region analysis ==="
REGION_LNC_DIR="$OUTDIR/lncseeker/per_region"
mkdir -p "$REGION_LNC_DIR"

# Load region info
REGION_GENES=()
REGION_STRINGS=()
while IFS=$'\t' read -r gene region; do
    REGION_GENES+=("$gene")
    REGION_STRINGS+=("$region")
done < <(python3 -c "
import json; c=json.load(open('$OUTDIR/config_resolved.json'))
for r in c.get('resolved_regions', []):
    print(f\"{r['gene']}\t{r['chrom']}:{r['start']}-{r['end']}\")
")

REGION_COUNT=${#REGION_GENES[@]}
if [ "$REGION_COUNT" -gt 0 ]; then
    for i in "${!REGION_GENES[@]}"; do
        gene="${REGION_GENES[$i]}"
        region="${REGION_STRINGS[$i]}"
        echo "  Exporting junctions for $gene ($region)"
        .venv/bin/python3 "$SCRIPT_DIR/export_lncseeker_junctions.py" "$CONFIG" \
            --region "$region" \
            -o "$REGION_LNC_DIR/${gene}_junctions.tsv"
    done
else
    echo "  (no gene regions configured)"
fi

echo ""
echo "=== Step 6: Per-region comparisons (genome-wide skipped) ==="
if [ "$REGION_COUNT" -gt 0 ]; then
    for i in "${!REGION_GENES[@]}"; do
        gene="${REGION_GENES[$i]}"
        region="${REGION_STRINGS[$i]}"
        region_dir="$OUTDIR/comparison/${gene}"
        mkdir -p "$region_dir"
        echo "  $gene ($region)"
        .venv/bin/python3 "$SCRIPT_DIR/compare_junctions.py" \
            --lncseeker "$REGION_LNC_DIR/${gene}_junctions.tsv" \
            --stringtie "$ST_JUNCTIONS_FILTERED" \
            --thresholds 1,2,5,10 \
            --region "$region" \
            -o "$region_dir/results.csv"
    done
else
    echo "  (no gene regions configured)"
fi

echo ""
echo "=== Step 8: Aggregated summary ==="
python3 -c "
import csv, json, sys, os

summary = []
comparison_dir = '$OUTDIR/comparison'
resolved = json.load(open('$OUTDIR/config_resolved.json'))
regions = resolved.get('resolved_regions', [])

for r in regions:
    gene = r['gene']
    csv_path = os.path.join(comparison_dir, gene, 'results.csv')
    if not os.path.exists(csv_path):
        continue
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            summary.append({
                'region': gene,
                'threshold': row['threshold'],
                'mode': row['mode'],
                'total_lncseeker': row['total_lncseeker'],
                'total_stringtie': row['total_stringtie'],
                'shared': row['shared'],
                'only_lncseeker': row['only_lncseeker'],
                'only_stringtie': row['only_stringtie'],
                'sensitivity': row['sensitivity'],
                'precision': row['precision'],
            })

# Append genome-wide
gw_path = os.path.join(comparison_dir, 'genome_wide', 'results.csv')
if os.path.exists(gw_path):
    with open(gw_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            summary.append({
                'region': 'genome_wide',
                'threshold': row['threshold'],
                'mode': row['mode'],
                'total_lncseeker': row['total_lncseeker'],
                'total_stringtie': row['total_stringtie'],
                'shared': row['shared'],
                'only_lncseeker': row['only_lncseeker'],
                'only_stringtie': row['only_stringtie'],
                'sensitivity': row['sensitivity'],
                'precision': row['precision'],
            })

out_path = os.path.join(comparison_dir, 'aggregated_summary.csv')
with open(out_path, 'w', newline='') as f:
    fields = ['region', 'threshold', 'mode', 'total_lncseeker', 'total_stringtie',
              'shared', 'only_lncseeker', 'only_stringtie', 'sensitivity', 'precision']
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(summary)
print(f'Aggregated summary: {out_path}', file=sys.stderr)
" 2>&1

echo ""
echo "=== Step 9: Classification against annotations (RefSeq + GENCODE lncRNA) ==="

REFSEQ_JUNCTIONS="$OUTDIR/stringtie/refseq_junctions.tsv"
GENCODE_JUNCTIONS="$OUTDIR/stringtie/gencode_junctions.tsv"

if [ ! -f "$REFSEQ_JUNCTIONS" ]; then
    echo "  Extracting RefSeq junctions..."
    .venv/bin/python3 "$SCRIPT_DIR/extract_gtf_junctions.py" "$GTF" -o "$REFSEQ_JUNCTIONS"
fi
GTF_GENCODE=$(python3 -c "import json; c=json.load(open('$OUTDIR/config_resolved.json')); print(c.get('gtf_gencode') or '')")
if [ -n "$GTF_GENCODE" ] && [ ! -f "$GENCODE_JUNCTIONS" ]; then
    echo "  Extracting GENCODE lncRNA junctions..."
    .venv/bin/python3 "$SCRIPT_DIR/extract_gtf_junctions.py" "$GTF_GENCODE" -o "$GENCODE_JUNCTIONS"
fi

CLASS_DIR="$OUTDIR/comparison/classification"
mkdir -p "$CLASS_DIR"

# Build common classification args
CLASS_ARGS=""
[ -f "$REFSEQ_JUNCTIONS" ] && CLASS_ARGS="$CLASS_ARGS --known-refseq $REFSEQ_JUNCTIONS"
[ -f "$GENCODE_JUNCTIONS" ] && CLASS_ARGS="$CLASS_ARGS --known-gencode $GENCODE_JUNCTIONS"

# Per-region classification (genome-wide skipped)
for i in "${!REGION_GENES[@]}"; do
    gene="${REGION_GENES[$i]}"
    region="${REGION_STRINGS[$i]}"
    .venv/bin/python3 "$SCRIPT_DIR/classify_junctions.py" \
        $CLASS_ARGS \
        --lncseeker "$REGION_LNC_DIR/${gene}_junctions.tsv" \
        --stringtie "$ST_JUNCTIONS" \
        --region "$region" \
        -o "$CLASS_DIR/${gene}.csv"
    echo ""
done

echo ""
echo "=== Step 10: Enriched aggregated summary ==="
python3 -c "
import csv, json, sys, os

comparison_dir = '$OUTDIR/comparison'
resolved = json.load(open('$OUTDIR/config_resolved.json'))
regions = resolved.get('resolved_regions', [])

# Load classification data
classification = {}
class_dir = os.path.join(comparison_dir, 'classification')
for fname in os.listdir(class_dir):
    if fname.endswith('.csv'):
        key = fname.replace('.csv', '')
        with open(os.path.join(class_dir, fname)) as f:
            reader = csv.DictReader(f)
            for row in reader:
                classification[key] = row

def enrich_row(row, cl):
    row['known_refseq_total'] = cl.get('known_refseq_total', '')
    row['known_gencode_total'] = cl.get('known_gencode_total', '')
    row['lnc_refseq_known'] = cl.get('lnc_refseq_known', '')
    row['lnc_gencode_known'] = cl.get('lnc_gencode_known', '')
    row['lnc_known_any'] = cl.get('lnc_known_any', '')
    row['lnc_novel'] = cl.get('lnc_novel', '')
    row['lnc_novel_also_st'] = cl.get('lnc_novel_also_st', '')
    row['lnc_novel_only'] = cl.get('lnc_novel_only', '')
    row['st_refseq_known'] = cl.get('st_refseq_known', '')
    row['st_gencode_known'] = cl.get('st_gencode_known', '')
    row['st_known_any'] = cl.get('st_known_any', '')
    row['st_recovered'] = cl.get('st_recovered', '')
    rtot = int(cl.get('known_refseq_total', 0) or 0)
    gtot = int(cl.get('known_gencode_total', 0) or 0)
    lnc_rs = int(cl.get('lnc_refseq_known', 0) or 0)
    lnc_gc = int(cl.get('lnc_gencode_known', 0) or 0)
    lnc_any = int(cl.get('lnc_known_any', 0) or 0)
    lnc_tot = int(row.get('total_lncseeker', 0) or 0)
    row['sens_lnc_vs_refseq'] = round(lnc_rs / rtot, 4) if rtot else ''
    row['prec_lnc_vs_refseq'] = round(lnc_rs / lnc_tot, 4) if lnc_tot else ''
    row['sens_lnc_vs_gencode'] = round(lnc_gc / gtot, 4) if gtot else ''
    row['prec_lnc_vs_gencode'] = round(lnc_gc / lnc_tot, 4) if lnc_tot else ''
    row['prec_lnc_vs_any'] = round(lnc_any / lnc_tot, 4) if lnc_tot else ''

summary = []
for r in regions:
    gene = r['gene']
    csv_path = os.path.join(comparison_dir, gene, 'results.csv')
    if not os.path.exists(csv_path):
        continue
    cl = classification.get(gene, {})
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            row['region'] = gene
            enrich_row(row, cl)
            summary.append(row)

out_path = os.path.join(comparison_dir, 'aggregated_summary.csv')
with open(out_path, 'w', newline='') as f:
    fields = ['region', 'threshold', 'mode',
              'total_lncseeker', 'total_stringtie', 'shared',
              'only_lncseeker', 'only_stringtie',
              'known_refseq_total', 'known_gencode_total',
              'lnc_refseq_known', 'lnc_gencode_known', 'lnc_known_any',
              'lnc_novel', 'lnc_novel_also_st', 'lnc_novel_only',
              'st_refseq_known', 'st_gencode_known', 'st_known_any',
              'st_recovered',
              'sensitivity', 'precision',
              'sens_lnc_vs_refseq', 'prec_lnc_vs_refseq',
              'sens_lnc_vs_gencode', 'prec_lnc_vs_gencode',
              'prec_lnc_vs_any']
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(summary)
print(f'Enriched aggregated summary: {out_path}', file=sys.stderr)
" 2>&1

echo ""
echo "=== Done ==="
echo "Results: $OUTDIR/comparison/"
echo "  <gene>/                 - per-region comparison"
echo "  classification/         - per-region classification"
echo "  aggregated_summary.csv  - enriched with novelty classification (RefSeq + GENCODE lncRNA)"
