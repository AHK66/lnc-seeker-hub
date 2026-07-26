#!/usr/bin/env python3
"""Compare splice junctions between lnc-seeker and StringTie.

Computes overlap, sensitivity, and precision at multiple read-count thresholds.
"""

import argparse
import csv
import sys


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome names to a common format (strip chr prefix, strip RefSeq NC_ prefix)."""
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


def load_junctions(path: str) -> list[tuple[str, int, int, str, int]]:
    """Load junctions from TSV: chrom, start, end, strand, count."""
    junctions = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = normalize_chrom(row['chrom'])
            start = int(row['start'])
            end = int(row['end'])
            strand = row['strand']
            count = int(row.get('reads', row.get('n_transcripts', 1)))
            junctions.append((chrom, start, end, strand, count))
    return junctions


def filter_by_region(
    junctions: list[tuple[str, int, int, str, int]],
    region_str: str,
) -> list[tuple[str, int, int, str, int]]:
    """Filter junctions to keep only those overlapping the given region.

    region_str format: chrom:start-end (0-based half-open).
    Junctions kept if they are fully within the region.
    """
    if not region_str:
        return junctions
    try:
        chrom_part, coord_part = region_str.split(':')
        rstart, rend = coord_part.split('-')
        rstart = int(rstart)
        rend = int(rend)
    except (ValueError, IndexError):
        print(f"Warning: could not parse region '{region_str}', ignoring filter", file=sys.stderr)
        return junctions

    chrom_part = normalize_chrom(chrom_part)

    result = []
    for chrom, start, end, strand, count in junctions:
        if chrom == chrom_part and start >= rstart and end <= rend:
            result.append((chrom, start, end, strand, count))
    return result


def build_sets(
    junctions: list[tuple[str, int, int, str, int]],
    threshold: int,
    stranded: bool,
) -> set:
    """Build set of junction keys filtered by count >= threshold."""
    result = set()
    for chrom, start, end, strand, count in junctions:
        if count < threshold:
            continue
        if stranded:
            result.add((chrom, start, end, strand))
        else:
            result.add((chrom, start, end))
    return result


def compare(
    lnc_path: str,
    st_path: str,
    thresholds: list[int],
    region: str = '',
) -> list[dict]:
    """Run comparison at each threshold (stranded and unstranded)."""
    lnc_juncs = load_junctions(lnc_path)
    st_juncs = load_junctions(st_path)

    if region:
        lnc_juncs = filter_by_region(lnc_juncs, region)
        st_juncs = filter_by_region(st_juncs, region)

    results = []

    for threshold in thresholds:
        for stranded in [True, False]:
            mode = 'stranded' if stranded else 'unstranded'

            lnc_set = build_sets(lnc_juncs, threshold, stranded)
            st_set = build_sets(st_juncs, threshold, stranded)

            shared = lnc_set & st_set
            only_lnc = lnc_set - st_set
            only_st = st_set - lnc_set

            sensitivity = len(shared) / len(st_set) if st_set else 0.0
            precision = len(shared) / len(lnc_set) if lnc_set else 0.0

            results.append({
                'threshold': threshold,
                'mode': mode,
                'total_lncseeker': len(lnc_set),
                'total_stringtie': len(st_set),
                'shared': len(shared),
                'only_lncseeker': len(only_lnc),
                'only_stringtie': len(only_st),
                'sensitivity': round(sensitivity, 4),
                'precision': round(precision, 4),
            })

    return results


def main():
    parser = argparse.ArgumentParser(description='Compare lnc-seeker vs StringTie junctions')
    parser.add_argument('--lncseeker', required=True, help='lnc-seeker junction TSV')
    parser.add_argument('--stringtie', required=True, help='StringTie junction TSV')
    parser.add_argument('-o', '--output', required=True, help='Output CSV path')
    parser.add_argument('--thresholds', default='1,2,5,10', help='Comma-separated thresholds')
    parser.add_argument('--region', default='', help='Filter to region (chrom:start-end, 0-based half-open)')
    args = parser.parse_args()

    thresholds = [int(t) for t in args.thresholds.split(',')]

    results = compare(args.lncseeker, args.stringtie, thresholds, region=args.region)

    fieldnames = [
        'threshold', 'mode', 'total_lncseeker', 'total_stringtie',
        'shared', 'only_lncseeker', 'only_stringtie', 'sensitivity', 'precision',
    ]

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    for r in results:
        print(
            f"thresh={r['threshold']:>2} {r['mode']:>10}: "
            f"lnc={r['total_lncseeker']:>5} st={r['total_stringtie']:>5} "
            f"shared={r['shared']:>5} only_lnc={r['only_lncseeker']:>5} "
            f"only_st={r['only_stringtie']:>5} "
            f"sens={r['sensitivity']:.4f} prec={r['precision']:.4f}"
        )


if __name__ == '__main__':
    main()
