#!/usr/bin/env python3
"""Classify junctions against RefSeq and GENCODE lncRNA annotations separately.

Two annotation sources (each optional):
  --known-refseq: RefSeq GTF junctions
  --known-gencode: GENCODE lncRNA GTF junctions

Output per source counts + combined novelty classification.
"""

import argparse
import csv
import os
import sys


def normalize_chrom(chrom: str) -> str:
    c = chrom
    if c.startswith('chr'):
        c = c[3:]
    if c.startswith('NC_'):
        parts = c.split('.')[0].split('_')
        if len(parts) >= 2:
            try:
                c = str(int(parts[-1]))
            except ValueError:
                pass
    return c


def load_junctions(path: str) -> set:
    result = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = normalize_chrom(row['chrom'])
            start = int(row['start'])
            end = int(row['end'])
            strand = row.get('strand', '.')
            result.add((chrom, start, end, strand))
    return result


def filter_by_region(junctions: set, region_str: str) -> set:
    if not region_str:
        return junctions
    try:
        chrom_part, coord_part = region_str.split(':')
        rstart, rend = coord_part.split('-')
        rstart, rend = int(rstart), int(rend)
    except (ValueError, IndexError):
        print(f"Warning: could not parse region '{region_str}'", file=sys.stderr)
        return junctions
    chrom_part = normalize_chrom(chrom_part)
    return {(c, s, e, strand) for (c, s, e, strand) in junctions
            if c == chrom_part and s >= rstart and e <= rend}


def classify(lnc_path: str, st_path: str,
             refseq_path: str | None, gencode_path: str | None,
             region: str = '') -> dict:
    lnc = load_junctions(lnc_path)
    st = load_junctions(st_path)

    refseq = load_junctions(refseq_path) if refseq_path else set()
    gencode = load_junctions(gencode_path) if gencode_path else set()

    known = refseq | gencode

    if region:
        known = filter_by_region(known, region)
        lnc = filter_by_region(lnc, region)
        st = filter_by_region(st, region)
        refseq = filter_by_region(refseq, region)
        gencode = filter_by_region(gencode, region)

    lnc_refseq_known = lnc & refseq
    lnc_gencode_known = lnc & gencode
    lnc_known_any = lnc & known
    lnc_novel = lnc - known
    lnc_novel_also_st = lnc_novel & (st - known)
    lnc_novel_only = lnc_novel - st

    st_refseq_known = st & refseq
    st_gencode_known = st & gencode
    st_known_any = st & known
    st_novel = st - known
    st_recovered = st_novel - lnc

    return {
        'region': region or 'genome_wide',

        'known_refseq_total': len(refseq),
        'known_gencode_total': len(gencode),

        'total_lncseeker': len(lnc),
        'lnc_refseq_known': len(lnc_refseq_known),
        'lnc_gencode_known': len(lnc_gencode_known),
        'lnc_known_any': len(lnc_known_any),
        'lnc_novel': len(lnc_novel),
        'lnc_novel_also_st': len(lnc_novel_also_st),
        'lnc_novel_only': len(lnc_novel_only),

        'total_stringtie': len(st),
        'st_refseq_known': len(st_refseq_known),
        'st_gencode_known': len(st_gencode_known),
        'st_known_any': len(st_known_any),
        'st_novel': len(st_novel),
        'st_recovered': len(st_recovered),
    }


def format_summary(s: dict) -> str:
    lines = [f"  RefSeq: {s['known_refseq_total']}, GENCODE lncRNA: {s['known_gencode_total']}"]
    lines.append(f"  lnc-seeker: {s['total_lncseeker']} total "
                 f"(RefSeq:{s['lnc_refseq_known']} Gencode:{s['lnc_gencode_known']} "
                 f"any:{s['lnc_known_any']}, novel:{s['lnc_novel']})")
    lines.append(f"    └ {s['lnc_novel_also_st']} also recovered by StringTie")
    lines.append(f"    └ {s['lnc_novel_only']} lnc-seeker only")
    lines.append(f"  StringTie: {s['total_stringtie']} total "
                 f"(RefSeq:{s['st_refseq_known']} Gencode:{s['st_gencode_known']} "
                 f"any:{s['st_known_any']}, novel:{s['st_novel']})")
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Classify junctions against RefSeq and GENCODE annotations')
    parser.add_argument('--known-refseq', help='RefSeq junctions TSV')
    parser.add_argument('--known-gencode', help='GENCODE lncRNA junctions TSV')
    parser.add_argument('--lncseeker', required=True, help='lnc-seeker junction TSV')
    parser.add_argument('--stringtie', required=True, help='StringTie junction TSV')
    parser.add_argument('-o', '--output', required=True, help='Output CSV path')
    parser.add_argument('--region', default='', help='Region filter (chrom:start-end)')
    args = parser.parse_args()

    if not args.known_refseq and not args.known_gencode:
        parser.error('At least one of --known-refseq or --known-gencode is required')

    stats = classify(args.lncseeker, args.stringtie,
                     args.known_refseq, args.known_gencode,
                     region=args.region)

    with open(args.output, 'w', newline='') as f:
        fields = ['region',
                  'known_refseq_total', 'known_gencode_total',
                  'total_lncseeker',
                  'lnc_refseq_known', 'lnc_gencode_known',
                  'lnc_known_any', 'lnc_novel',
                  'lnc_novel_also_st', 'lnc_novel_only',
                  'total_stringtie',
                  'st_refseq_known', 'st_gencode_known',
                  'st_known_any', 'st_novel', 'st_recovered']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(stats)

    print(format_summary(stats), file=sys.stderr)


if __name__ == '__main__':
    main()
