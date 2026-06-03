#!/usr/bin/env python3
"""
merge_counts.py

Merge per-sample featureCounts output files into a single count matrix
suitable for DESeq2 input.

featureCounts output format (tab-separated):
  # Comment line(s)
  Geneid  Chr  Start  End  Strand  Length  /path/to/sample.bam
  gene1   1    100    200  +       100     42
  ...

Usage:
  merge_counts.py --input s1.featureCounts.txt s2.featureCounts.txt ... --output merged_counts.txt
  merge_counts.py --input *.featureCounts.txt --output merged_counts.txt
"""

import argparse
import os
import sys
import glob


def parse_featurecounts(filepath):
    """
    Parse one featureCounts .txt file.
    Returns (gene_ids list, sample_name string, counts list).
    """
    gene_ids    = []
    counts      = []
    sample_name = None

    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Skip comment / Program line
            if line.startswith('#'):
                continue

            cols = line.split('\t')

            # Header row — last column is the BAM path
            if cols[0] == 'Geneid':
                bam_path = cols[-1]
                # Strip directory and extension, clean up STAR suffixes
                base = os.path.basename(bam_path)
                base = os.path.splitext(base)[0]
                for suffix in ['.Aligned.sortedByCoord.out', '.sorted', '.bam']:
                    if base.endswith(suffix):
                        base = base[: -len(suffix)]
                sample_name = base
                continue

            # Data rows — col 0 = gene_id, last col = count
            gene_ids.append(cols[0])
            counts.append(cols[-1])   # keep as string; DESeq2 needs integers

    if sample_name is None:
        raise ValueError(f"Could not find header row in {filepath}")

    return gene_ids, sample_name, counts


def main():
    parser = argparse.ArgumentParser(
        description='Merge featureCounts files into a single count matrix.'
    )
    parser.add_argument(
        '--input', nargs='+', required=True,
        help='featureCounts .txt files (space-separated or glob)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output merged count matrix (tab-separated)'
    )
    args = parser.parse_args()

    # Expand any globs that the shell didn't expand (e.g. when called from Nextflow)
    input_files = []
    for pattern in args.input:
        expanded = glob.glob(pattern)
        if expanded:
            input_files.extend(expanded)
        else:
            input_files.append(pattern)   # let open() raise a clear error

    # Keep only .txt files and sort for reproducibility
    input_files = sorted([f for f in input_files if f.endswith('.txt')])

    if not input_files:
        sys.exit('ERROR: No .txt featureCounts files found in --input.')

    print(f'Merging {len(input_files)} featureCounts file(s)...')

    all_data   = {}     # {sample_name: [count_strings]}
    gene_order = None

    for fpath in input_files:
        if not os.path.isfile(fpath):
            sys.exit(f'ERROR: File not found: {fpath}')

        gene_ids, sample_name, counts = parse_featurecounts(fpath)

        if gene_order is None:
            gene_order = gene_ids
        elif gene_order != gene_ids:
            sys.exit(
                f'ERROR: Gene order / count in {fpath} does not match the '
                f'first file. All featureCounts files must use the same GTF.'
            )

        if sample_name in all_data:
            sys.exit(f'ERROR: Duplicate sample name "{sample_name}". '
                     'Check that all BAM files have unique base names.')

        all_data[sample_name] = counts

    samples = sorted(all_data.keys())

    with open(args.output, 'w') as out:
        # Header
        out.write('gene_id\t' + '\t'.join(samples) + '\n')
        # Data rows
        for i, gene in enumerate(gene_order):
            row = [gene] + [all_data[s][i] for s in samples]
            out.write('\t'.join(row) + '\n')

    print(
        f'Done. Matrix: {len(gene_order)} genes × {len(samples)} samples'
        f'  →  {args.output}'
    )


if __name__ == '__main__':
    main()
