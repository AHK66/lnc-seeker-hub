#!/usr/bin/env python3
"""Extract known exon-exon junctions from a GTF annotation file.

Output: TSV with columns chrom, start, end, strand, transcript_id, gene_name.

Junctions are 0-based half-open (same as StringTie and lnc-seeker).
A junction between exon_i and exon_{i+1} of the same transcript:
  start = exon_i.end    (1-based GTF exon end → 0-based intron start)
  end   = exon_{i+1}.start - 1  (1-based GTF exon start → 0-based intron end)
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict


ATTR_RE = re.compile(r'(gene|transcript_id|gene_name)\s+"([^"]+)"')


def parse_attributes(attr_str: str) -> dict:
    return dict(ATTR_RE.findall(attr_str))


def extract_junctions(gtf_path: str) -> list[dict]:
    exons_by_transcript = defaultdict(list)
    open_fn = gzip.open if gtf_path.endswith('.gz') else open

    with open_fn(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'exon':
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            transcript_id = attrs.get('transcript_id', '')
            if not transcript_id:
                continue
            exon_num = int(attrs.get('exon_number', 0))
            exons_by_transcript[transcript_id].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'exon_number': exon_num,
                'gene': attrs.get('gene_name', attrs.get('gene', '')),
            })

    junctions = []
    for tid, exons in exons_by_transcript.items():
        exons.sort(key=lambda e: e['start'])
        for i in range(len(exons) - 1):
            e1 = exons[i]
            e2 = exons[i + 1]
            if e1['chrom'] != e2['chrom'] or e1['strand'] != e2['strand']:
                continue
            jstart = e1['end']
            jend = e2['start'] - 1
            if jstart >= jend:
                continue
            junctions.append({
                'chrom': e1['chrom'],
                'start': jstart,
                'end': jend,
                'strand': e1['strand'],
                'transcript_id': tid,
                'gene': e1['gene'],
            })

    return junctions


def main():
    parser = argparse.ArgumentParser(description='Extract exon-exon junctions from GTF(s)')
    parser.add_argument('gtf', nargs='+', help='GTF file(s) (can be .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV path')
    args = parser.parse_args()

    all_junctions = []
    seen = set()
    for gtf_path in args.gtf:
        junctions = extract_junctions(gtf_path)
        for j in junctions:
            key = (j['chrom'], j['start'], j['end'], j['strand'])
            if key not in seen:
                seen.add(key)
                all_junctions.append(j)
        print(f'  {len(junctions)} from {gtf_path}', file=sys.stderr)

    all_junctions.sort(key=lambda j: (j['chrom'], j['start'], j['end']))

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['chrom', 'start', 'end', 'strand', 'transcript_id', 'gene'])
        for j in all_junctions:
            writer.writerow([j['chrom'], j['start'], j['end'], j['strand'],
                             j['transcript_id'], j['gene']])

    print(f'Extracted {len(all_junctions)} unique known junctions from {len(args.gtf)} GTF files', file=sys.stderr)


if __name__ == '__main__':
    main()
