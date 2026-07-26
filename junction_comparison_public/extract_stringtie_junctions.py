#!/usr/bin/env python3
"""Extract splice junctions from a StringTie GTF output file.

StringTie GTF has exon features grouped by transcript_id.
Junctions are the gaps between consecutive exons within each transcript.
StringTie coordinates are 1-based; output is converted to 0-based half-open.
"""

import argparse
import csv
import sys
from collections import defaultdict


def parse_attr(attr_str: str) -> dict:
    """Parse GTF attribute column into a dict."""
    result = {}
    i = 0
    while i < len(attr_str):
        while i < len(attr_str) and attr_str[i] in (' ', '\t', ';'):
            i += 1
        if i >= len(attr_str):
            break
        end = attr_str.index(' ', i)
        key = attr_str[i:end]
        i = end + 1
        if i < len(attr_str) and attr_str[i] == '"':
            i += 1
            end = attr_str.index('"', i)
            value = attr_str[i:end]
            i = end + 2
        elif i < len(attr_str) and attr_str[i] != ';':
            end = attr_str.index(';', i) if ';' in attr_str[i:] else len(attr_str)
            value = attr_str[i:end].strip()
            i = end + 1
        else:
            value = ''
            i += 1
        result[key] = value
    return result


def extract_junctions(gtf_path: str) -> list[tuple[str, int, int, str, int]]:
    """Extract junctions from StringTie GTF.

    Returns list of (chrom, start_0based, end_0based, strand, n_exon_pairs).
    """
    exons_by_transcript: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom = parts[0]
            feature = parts[2]
            if feature != 'exon':
                continue
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attr(parts[8])
            tid = attrs.get('transcript_id', '')
            exons_by_transcript[(chrom, tid)].append((start, end, strand))

    junction_set: dict[tuple[str, int, int, str], int] = defaultdict(int)

    for (chrom, tid), exon_list in exons_by_transcript.items():
        exon_list.sort(key=lambda x: x[0])
        strand = exon_list[0][2] if exon_list else '.'
        for i in range(len(exon_list) - 1):
            curr_end = exon_list[i][1]
            next_start = exon_list[i + 1][0]
            gap = next_start - curr_end - 1
            if gap <= 0:
                continue
            jstart = curr_end
            jend = next_start - 1
            key = (chrom, jstart, jend, strand)
            junction_set[key] += 1

    result = [(chrom, s, e, strand, count) for (chrom, s, e, strand), count in sorted(junction_set.items())]
    return result


def main():
    parser = argparse.ArgumentParser(description='Extract junctions from StringTie GTF')
    parser.add_argument('gtf', help='StringTie output GTF')
    parser.add_argument('-o', '--output', required=True, help='Output TSV path')
    args = parser.parse_args()

    junctions = extract_junctions(args.gtf)

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['chrom', 'start', 'end', 'strand', 'n_transcripts'])
        for chrom, start, end, strand, count in junctions:
            writer.writerow([chrom, start, end, strand, count])

    print(f'Extracted {len(junctions)} unique junctions from StringTie GTF', file=sys.stderr)


if __name__ == '__main__':
    main()
