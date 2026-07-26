#!/usr/bin/env python3
"""Run lnc-seeker analysis and export aggregated cohort-level junctions to TSV.

Accepts junction_comparison config.cfg format directly and generates the
lnc-seeker JSON config internally.
"""

import argparse
import csv
import json
import os
import sys

# Ensure we can import from junction_comparison
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lnc_seeker
from parse_config import parse_config, resolve_gene_regions


def aggregate_junctions(samples: dict) -> list[dict]:
    """Aggregate junctions across all samples in a cohort.

    For each unique (chrom, start, end, strand), sum reads and reads_clean.
    """
    from collections import defaultdict
    agg: dict[tuple[str, int, int, str], dict] = defaultdict(lambda: {'reads': 0, 'reads_clean': 0})

    for sample_name, sample_data in samples.items():
        spans = sample_data.get('junction_spans', []) if isinstance(sample_data, dict) else []
        for span in spans:
            key = (span['reference'], span['start'], span['end'], span['junction_strand'])
            agg[key]['reads'] += span['reads']
            agg[key]['reads_clean'] += span.get('reads_clean', 0)

    result = []
    for (chrom, start, end, strand), counts in sorted(agg.items()):
        result.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'reads': counts['reads'],
            'reads_clean': counts['reads_clean'],
        })
    return result


def build_analysis_config(parsed_config: dict, region: str | None = None) -> dict:
    """Build lnc-seeker analysis JSON config from parsed comparison config."""
    if region:
        chrom, coords = region.split(':', 1)
        start_s, end_s = coords.split('-', 1)
        analysis_reference = chrom
        analysis_start = int(start_s)
        analysis_end = int(end_s)
    else:
        analysis_reference = None
        analysis_start = None
        analysis_end = None

    gtf_paths = [parsed_config['gtf']]
    gencode = parsed_config.get('gtf_gencode')
    if gencode:
        gtf_paths.append(gencode)

    config = {
        'data_selection': {
            'bam_paths': parsed_config['bam_files'],
            'gtf_paths': gtf_paths,
            'genome_path': parsed_config['genome'],
            'filter_annotations': False,
            'max_cache_memory_mb': 65536.0,
        },
        'coverage_and_junctions_profile': {
            'min_mapping_quality': 20,
        },
    }

    if analysis_reference is not None:
        config['data_selection']['analysis_reference'] = analysis_reference
        config['data_selection']['analysis_start'] = analysis_start
        config['data_selection']['analysis_end'] = analysis_end

    if parsed_config.get('assembly_report_path'):
        config['data_selection']['assembly_report_path'] = parsed_config['assembly_report_path']

    if parsed_config.get('cache_dir'):
        config['data_selection']['cache_dir'] = parsed_config['cache_dir']

    return config


def main():
    parser = argparse.ArgumentParser(description='Export aggregated lnc-seeker junctions')
    parser.add_argument('config', help='junction_comparison config.cfg path')
    parser.add_argument('-o', '--output', help='Output aggregated junction TSV path (default: <output_dir>/lncseeker_junctions.tsv)')
    parser.add_argument('--region', help='Restrict analysis to a genomic region (chrom:start-end)')
    args = parser.parse_args()

    parsed = parse_config(args.config)

    output_path = args.output or os.path.join(parsed['output_dir'], 'lncseeker_junctions.tsv')
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    if not parsed['bam_files']:
        print('No BAM files found in config', file=sys.stderr)
        return 1

    analysis_config = build_analysis_config(parsed, region=args.region)

    progress = lnc_seeker.SessionProgress()
    result_json = lnc_seeker.run_analysis_py(json.dumps(analysis_config), progress)
    result = json.loads(result_json)

    samples = result.get('samples', {})
    if not samples:
        print('No samples found in analysis result', file=sys.stderr)
        return 1

    junctions = aggregate_junctions(samples)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['chrom', 'start', 'end', 'strand', 'reads', 'reads_clean'])
        for j in junctions:
            writer.writerow([j['chrom'], j['start'], j['end'], j['strand'], j['reads'], j['reads_clean']])

    print(f'Aggregated {len(junctions)} unique junctions from {len(samples)} samples to {output_path}', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
