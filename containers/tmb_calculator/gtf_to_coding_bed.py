#!/usr/bin/env python3
"""
Extract coding regions from a GENCODE GTF (e.g., gencode.v46.annotation.gtf.gz)
and write a merged BED of CDS intervals for protein-coding transcripts, per gene.

Output is 0-based, half-open BED intervals (chrom, start, end, gene_name), merged per gene.

Examples:

  python gtf_to_coding_bed.py gencode.v46.annotation.gtf.gz -o coding_regions.bed
  python gtf_to_coding_bed.py gencode.v46.annotation.gtf.gz -o coding_regions.bed.gz

Notes:
- Uses feature == "CDS" and transcript_type (or transcript_biotype) == "protein_coding".
- Merges CDS per gene, outputs gene name in 4th column.
- Limits to primary chromosomes chr1..chr22, chrX, chrY by default (excludes chrM and alt contigs).
  Use --include-chrM to include chrM.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterable

import click


PrimaryChroms = {f"chr{i}" for i in range(1, 23)}.union({"chrX", "chrY"})


def open_maybe_gzip(path: str):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    """Parse the attributes column of a GTF line into a dict.

    Handles formats like: key "value"; key2 "value2";
    Also tolerates key=value; style if present.
    """
    attrs: Dict[str, str] = {}
    # First try key "value"; patterns
    for m in re.finditer(r'(\S+)\s+"(.*?)"\s*;?', attr):
        attrs[m.group(1)] = m.group(2)
    # Also pick up key=value (no quotes) patterns
    for m in re.finditer(r'(\S+?)=(.*?)(?:;|$)', attr):
        key = m.group(1)
        val = m.group(2).strip().strip('"')
        if key not in attrs:
            attrs[key] = val
    return attrs


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged: List[Tuple[int, int]] = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged


@click.command()
@click.argument('gtf', type=click.Path(exists=True, dir_okay=False))
@click.option('-o', '--output', required=True, help='Output BED file (gz accepted)')
@click.option('--include-chrM', is_flag=True, default=False, help='Include chrM (mitochondrial) in output')
def main(gtf: str, output: str, include_chrm: bool):
    # Collect intervals per chromosome and gene
    intervals_by_chr_gene: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}

    include_set = set(PrimaryChroms)
    if include_chrm:
        include_set.add('chrM')

    with open_maybe_gzip(gtf) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr = parts

            # Keep only primary chromosomes (exclude alts/patches by default)
            if chrom not in include_set:
                continue

            # Only coding sequence features
            if feature != 'CDS':
                continue

            try:
                start_1 = int(start)
                end_1 = int(end)
            except ValueError:
                continue
            if end_1 < start_1:
                continue

            attrs = parse_gtf_attributes(attr)
            ttype = attrs.get('transcript_type') or attrs.get('transcript_biotype')
            if ttype != 'protein_coding':
                continue

            gene_name = (
                attrs.get('gene_name')
                or attrs.get('gene')
                or attrs.get('gene_id')
                or 'NA'
            )

            # Convert to 0-based half-open
            s0 = start_1 - 1
            e0 = end_1
            intervals_by_chr_gene.setdefault(chrom, {}).setdefault(gene_name, []).append((s0, e0))

    # Merge and write
    out_path = Path(output)
    opener = gzip.open if out_path.suffix == '.gz' else open
    mode = 'wt'

    with opener(out_path, mode) as out:
        for chrom in sorted(intervals_by_chr_gene.keys(), key=lambda c: (c not in PrimaryChroms, c)):
            genes = intervals_by_chr_gene[chrom]
            for gene in sorted(genes.keys()):
                merged = merge_intervals(genes[gene])
                for s, e in merged:
                    out.write(f"{chrom}\t{s}\t{e}\t{gene}\n")

    click.echo(f"Wrote coding BED (with gene names) to {output}")


if __name__ == '__main__':
    main()
