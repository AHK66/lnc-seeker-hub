#!/usr/bin/env python3
"""Parse the dedicated junction comparison config.cfg file.

Format:
    stringtie=/usr/local/bin/stringtie
    gtf=/path/to/annotation.gtf.gz
    genome=/path/to/genome.fa
    gene_region=MALAT1:500
    bam=/path/to/sample1.bam
    bam=/path/to/sample2.bam
"""

import argparse
import gzip
import json
import sys


def parse_config(path: str) -> dict:
    result = {
        'stringtie': None,
        'gtf': None,
        'gtf_stringtie': None,
        'gtf_gencode': None,
        'genome': None,
        'output_dir': 'junction_comparison_output',
        'cache_dir': None,
        'gene_regions': [],
        'bam_files': [],
        'resolved_regions': [],
    }
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.split('#')[0].strip()
            if key == 'stringtie':
                result['stringtie'] = value
            elif key == 'gtf':
                result['gtf'] = value
            elif key == 'gtf_stringtie':
                result['gtf_stringtie'] = value
            elif key == 'gtf_gencode':
                result['gtf_gencode'] = value
            elif key == 'genome':
                result['genome'] = value
            elif key == 'output_dir':
                result['output_dir'] = value
            elif key == 'cache_dir':
                result['cache_dir'] = value
            elif key == 'gene_region':
                result['gene_regions'].append(value)
            elif key == 'bam':
                result['bam_files'].append(value)
    return result


def resolve_gene_regions(gtf_path: str, gene_region_specs: list[str], padding_default: int = 5000) -> list[dict]:
    regions = []
    specs = []
    for spec in gene_region_specs:
        if ':' in spec:
            gene_name, padding_str = spec.rsplit(':', 1)
            try:
                padding = int(padding_str)
            except ValueError:
                padding = padding_default
        else:
            gene_name = spec
            padding = padding_default
        specs.append((gene_name, padding))

    if not specs or not gtf_path:
        return regions

    open_fn = gzip.open if gtf_path.endswith('.gz') else open
    try:
        with open_fn(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'gene':
                    continue
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                attrs = parts[8]
                for gene_name, padding in specs:
                    if (f'gene_name "{gene_name}"' in attrs or f'gene_name "{gene_name};' in attrs
                        or f'gene "{gene_name}"' in attrs or f'gene "{gene_name};' in attrs):
                        regions.append({
                            'gene': gene_name,
                            'chrom': chrom,
                            'start': max(0, start - padding),
                            'end': end + padding,
                            'padding': padding,
                        })
                        break
    except FileNotFoundError:
        print(f"Warning: GTF not found at {gtf_path}", file=sys.stderr)

    return regions


def main():
    parser = argparse.ArgumentParser(description='Parse junction comparison config')
    parser.add_argument('config', help='Path to config_comparison.cfg')
    parser.add_argument('-o', '--output', help='Output JSON path (optional)')
    parser.add_argument('--resolve-genes', action='store_true',
                        help='Resolve gene_region names to coordinates from GTF')
    args = parser.parse_args()

    config = parse_config(args.config)

    if args.resolve_genes and config['gtf']:
        config['resolved_regions'] = resolve_gene_regions(config['gtf'], config['gene_regions'])

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(config, f, indent=2)

    print(f'StringTie: {config["stringtie"]}', file=sys.stderr)
    print(f'GTF (lnc-seeker): {config["gtf"]}', file=sys.stderr)
    print(f'GTF (StringTie): {config.get("gtf_stringtie", config["gtf"])}', file=sys.stderr)
    print(f'Genome: {config["genome"]}', file=sys.stderr)
    print(f'Gene region specs: {len(config["gene_regions"])}', file=sys.stderr)
    print(f'Resolved regions: {len(config["resolved_regions"])}', file=sys.stderr)
    print(f'BAM files: {len(config["bam_files"])}', file=sys.stderr)
    print(json.dumps(config, indent=2))


if __name__ == '__main__':
    main()
