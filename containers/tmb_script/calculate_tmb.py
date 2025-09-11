#!/usr/bin/env python3
"""
TMB Calculator for VEP-annotated VCF files
Calculates Tumor Mutational Burden
"""

import gzip
import json
import sys
import re
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any
from bisect import bisect_left

import click
from cyvcf2 import VCF
from rich.console import Console
from rich.progress import track
from rich.table import Table

console = Console()

# Nonsynonymous consequence types
NONSYNONYMOUS_CONSEQUENCES = {
    'missense_variant',
    'frameshift_variant',
    'stop_gained',
    'stop_lost',
    'start_lost',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'inframe_insertion',
    'inframe_deletion',
    'protein_altering_variant',
    'feature_elongation',
    'feature_truncation',
    'incomplete_terminal_codon_variant',
    'transcript_truncation'
}

SYNONYMOUS_CONSEQUENCES = {'synonymous_variant', 'stop_retained_variant', 'start_retained_variant'}

# All coding consequence types (for TMB calculation)
CODING_CONSEQUENCES = NONSYNONYMOUS_CONSEQUENCES.union(SYNONYMOUS_CONSEQUENCES)
class VEPParser:
    """Parses VEP CSQ header and annotation strings (including MAX_AF)."""

    def __init__(self, fields: List[str]):
        self.fields = fields
        self.index = {name: i for i, name in enumerate(fields)}

        # Validate required fields
        if 'Consequence' not in self.index:
            raise RuntimeError("VEP CSQ header missing required 'Consequence' field")
        if 'MAX_AF' not in self.index:
            # Per requirement: must be present, otherwise fail
            raise RuntimeError("VEP CSQ header missing required 'MAX_AF' field from gnomAD annotations")

    @classmethod
    def from_vcf_header(cls, raw_header: str) -> 'VEPParser':
        """Extract CSQ Format fields from VCF raw header and return a parser instance."""
        csq_lines = [line for line in raw_header.splitlines() if line.startswith('##INFO=<ID=CSQ')]
        if not csq_lines:
            raise RuntimeError("VCF is missing VEP CSQ INFO header")

        line = csq_lines[0]
        # Capture the 'Format: ...' content
        m = re.search(r'Format: (.+?)">', line)
        if m:
            fmt = m.group(1)
        else:
            raise RuntimeError("Unable to parse VEP CSQ header format")
        fields = fmt.strip().split('|')
        return cls(fields)

    def parse(self, csq_string: str) -> List[List[str]]:
        if not csq_string:
            return []
        # This is used in case of multiple transcripts for the same variant
        # but we don't actually use it because the input is assumed to be single-transcript
        annotations = csq_string.split(',')
        parsed: List[List[str]] = []
        for ann in annotations:
            parts = ann.split('|')
            # This should not happen, as each annotation should have the same number of fields as the header format
            if len(parts) < len(self.fields):
                raise RuntimeError(f"Invalid CSQ length after splitting: {csq_string} with fields: {self.fields}")
            parsed.append(parts)
        return parsed

    def get_consequences(self, csq_string: str) -> List[str]:
        results: List[str] = []
        idx = self.index['Consequence']
        for parts in self.parse(csq_string):
            consequence = parts[idx] if len(parts) > idx else ''
            if consequence:
                results.extend([c for c in consequence.split('&') if c])
        return results

    def get_max_af(self, csq_string: str) -> float:
        idx = self.index['MAX_AF']
        max_val = 0.0
        seen = False
        for parts in self.parse(csq_string):
            val_str = parts[idx] if len(parts) > idx else ''
            if val_str in (None, '', '.'):
                val = 0.0
            else:
                try:
                    val = float(val_str)
                except ValueError:
                    val = 0.0
            seen = True
            if val > max_val:
                max_val = val
        return max_val if seen else 0.0


class TMBCalculator:
    def __init__(self, vcf_path: str, coverage_path: str, min_vaf: float = 0.10,
                 min_depth: int = 10, min_coverage: int = 10, show_consequences: bool = False,
                 max_af_threshold: float = 0.03, debug_tsv: Optional[str] = None,
                 region_bed: Optional[str] = None):
        self.vcf_path = vcf_path
        self.coverage_path = coverage_path
        self.min_vaf = min_vaf
        self.min_depth = min_depth
        self.min_coverage = min_coverage
        self.show_consequences = show_consequences
        self.max_af_threshold = max_af_threshold  # gnomAD MAX_AF threshold
        self.vep_parser: Optional[VEPParser] = None
        self.debug_tsv = debug_tsv
        self._debug_rows: List[Dict[str, Any]] = []
        self.region_bed_path = region_bed
        # Interval index built from BED: chrom -> list of (start, end) merged, sorted
        self._bed_intervals: Dict[str, List[Tuple[int, int]]] = {}
        # For binary search acceleration: chrom -> list of starts
        self._bed_starts: Dict[str, List[int]] = {}

        self.stats = {
            'total_variants': 0,
            'pass_variants': 0,
            'coding_variants': 0,
            'nonsynonymous_variants': 0,
            'synonymous_variants': 0,
            'noncoding_variants': 0,
            'filtered_all_variants': 0,
            'filtered_coding_variants': 0,
            'filtered_nonsynonymous_variants': 0,
            'eligible_bases': 0,
            'eligible_mbp': 0,
            'tmb_score': 0,
            'nonsynonymous_score': 0,
            # Denominator used for nonsynonymous score (Mbp)
            'nonsynonymous_denominator_mbp': 0.0,
            # If no region BED is provided, we assume whole-genome VCF and CDS size of 34 Mbp
            'assumed_cds_mbp': None,
            'assumed_vcf_whole_genome': False,
            'average_coverage': 0,
            'coding_consequence_counts': {},
            'noncoding_consequence_counts': {},
            'filtered_coding_consequence_counts': {},
            'filtered_noncoding_consequence_counts': {},
            # Region BED metadata
            'region_bed_used': False,
            'region_bed_path': '',
            'region_bed_total_bases': 0,
            'region_bed_total_mbp': 0.0,
        }

        if self.region_bed_path:
            self._load_region_bed(self.region_bed_path)

    # -------------------------
    # Region BED handling
    # -------------------------
    def _load_region_bed(self, bed_path: str):
        console.print(f"[blue]Loading region BED from {bed_path}...")
        if not Path(bed_path).exists():
            raise FileNotFoundError(f"Region BED not found: {bed_path}")
        open_func = gzip.open if bed_path.endswith('.gz') else open
        intervals_by_chr: Dict[str, List[Tuple[int, int]]] = {}
        with open_func(bed_path, 'rt') as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    console.print(f"[yellow]Warning: Line {line.strip()} in BED file has fewer than 3 fields, skipping")
                    continue
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                # Skip mitochondrial by default (consistent with rest of tool)
                if chrom in ['chrM', 'MT', 'chrMT', 'M']:
                    continue
                intervals_by_chr.setdefault(chrom, []).append((start, end))

        # Merge and index
        def merge(iv: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
            if not iv:
                return []
            iv.sort(key=lambda t: t[0])
            merged: List[Tuple[int, int]] = []
            cs, ce = iv[0]
            for s, e in iv[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    merged.append((cs, ce))
                    cs, ce = s, e
            merged.append((cs, ce))
            return merged

        for chrom, iv in intervals_by_chr.items():
            merged = merge(iv)
            self._bed_intervals[chrom] = merged
            self._bed_starts[chrom] = [s for s, _ in merged]

        total_bed_bases = sum(e - s for iv in self._bed_intervals.values() for s, e in iv)
        self.stats['region_bed_used'] = True
        self.stats['region_bed_path'] = bed_path
        self.stats['region_bed_total_bases'] = total_bed_bases
        self.stats['region_bed_total_mbp'] = total_bed_bases / 1_000_000
        console.print(f"[green]Loaded {len(self._bed_intervals)} chromosomes, {total_bed_bases:,} total bases in BED")

    def _overlap_length(self, chrom: str, start: int, end: int) -> int:
        """Return overlap length between [start, end) and BED intervals for chrom.
        If no BED provided, returns full length.
        """
        if not self._bed_intervals:
            return max(0, end - start)
        iv = self._bed_intervals.get(chrom)
        if not iv:
            return 0
        starts = self._bed_starts[chrom]
        # Find first interval with start <= end via binary search
        i = bisect_left(starts, start)
        # Step back one to catch interval that starts before but overlaps
        if i > 0:
            i -= 1
        overlap = 0
        n = len(iv)
        while i < n and iv[i][0] < end:
            s, e = iv[i]
            if e <= start:
                i += 1
                continue
            # overlap with current interval
            os = max(s, start)
            oe = min(e, end)
            if oe > os:
                overlap += (oe - os)
            if e >= end:
                break
            i += 1
        return overlap

    def _pos_in_bed(self, chrom: str, pos0: int) -> bool:
        """Check if 0-based position lies within any BED interval for chrom."""
        if not self._bed_intervals:
            return True
        iv = self._bed_intervals.get(chrom)
        if not iv:
            return False
        starts = self._bed_starts[chrom]
        i = bisect_left(starts, pos0)
        if i == len(iv):
            i -= 1
        # Check current and previous interval for containment
        for j in (i, i - 1):
            if 0 <= j < len(iv):
                s, e = iv[j]
                if s <= pos0 < e:
                    return True
        return False

    def parse_vep_consequences(self, csq_string: str) -> List[str]:
        """Parse VEP CSQ field using header-defined indices and extract consequence types"""
        if not self.vep_parser:
            raise RuntimeError("VEP parser not initialized before parsing consequences")
        return self.vep_parser.get_consequences(csq_string)

    def is_nonsynonymous(self, consequences: List[str]) -> bool:
        """Check if variant has nonsynonymous consequences"""
        return any(cons in NONSYNONYMOUS_CONSEQUENCES for cons in consequences)

    def is_coding(self, consequences: List[str]) -> bool:
        """Check if variant has coding consequences (synonymous or nonsynonymous)"""
        return any(cons in CODING_CONSEQUENCES for cons in consequences)

    def is_synonymous(self, consequences: List[str]) -> bool:
        """Check if variant has synonymous consequences"""
        return any(cons in SYNONYMOUS_CONSEQUENCES for cons in consequences)

    def calculate_eligible_regions(self) -> Tuple[int, float]:
        """Calculate eligible regions from mosdepth coverage file"""
        console.print(f"[blue]Analyzing coverage from {self.coverage_path}...")

        eligible_bases = 0
        total_bases = 0
        total_coverage = 0
        region_count = 0

        open_func = gzip.open if self.coverage_path.endswith('.gz') else open

        with open_func(self.coverage_path, 'rt') as f:
            for line in track(f, description="Processing coverage regions"):
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    console.print(f"[yellow]Warning: Line {line.strip()} in coverage file has fewer than 4 fields, skipping")
                    continue

                try:
                    chrom, start, end, coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])

                    # Skip mitochondrial chromosome
                    if chrom in ['chrM', 'MT', 'chrMT', 'M']:
                        continue

                    # If a region BED is provided, only consider the overlap
                    if self._bed_intervals:
                        overlap_len = self._overlap_length(chrom, start, end)
                        if overlap_len <= 0:
                            continue
                        total_bases += overlap_len
                        total_coverage += coverage * overlap_len
                        region_count += 1
                        if coverage >= self.min_coverage:
                            eligible_bases += overlap_len
                    else:
                        region_size = end - start
                        total_bases += region_size
                        total_coverage += coverage * region_size
                        region_count += 1
                        if coverage >= self.min_coverage:
                            eligible_bases += region_size

                except (ValueError, IndexError):
                    continue

        average_coverage = total_coverage / total_bases if total_bases > 0 else 0
        self.stats['average_coverage'] = average_coverage
        self.stats['eligible_bases'] = eligible_bases
        self.stats['eligible_mbp'] = eligible_bases / 1_000_000

        console.print(f"[green]Processed {region_count:,} coverage regions")
        console.print(f"[green]Average coverage: {average_coverage:.2f}x")
        console.print(f"[green]Eligible bases (≥{self.min_coverage}x): {eligible_bases:,} ({eligible_bases/1_000_000:.2f} Mbp)")

        return eligible_bases, average_coverage

    def process_variants(self) -> Tuple[int, int, int]:
        """Process VCF file and count qualifying variants"""
        console.print(f"[blue]Processing variants from {self.vcf_path}...")

        vcf = VCF(self.vcf_path)
        # Enforce single-sample VCF
        samples = getattr(vcf, 'samples', [])
        if len(samples) != 1:
            console.print(f"[red]Error: VCF contains {len(samples)} samples; expected exactly 1.")
            sys.exit(1)
        # Initialize VEP parser from VCF header (requires MAX_AF)
        self.vep_parser = VEPParser.from_vcf_header(vcf.raw_header)
        qualifying_all_variants = 0
        qualifying_coding_variants = 0
        qualifying_nonsynonymous_variants = 0

        for variant in track(vcf, description="Processing variants"):
            # Skip mitochondrial variants (applies to all counts)
            if variant.CHROM in ['chrM', 'MT', 'chrMT', 'M']:
                continue

            # Restrict to region BED if provided (applies to all counts)
            if self._bed_intervals and not self._pos_in_bed(variant.CHROM, variant.POS - 1):
                continue

            # Count total variants after region filtering
            self.stats['total_variants'] += 1

            # Now enforce PASS for downstream consequence/calling stats
            if variant.FILTER is not None and variant.FILTER != 'PASS':
                continue
            self.stats['pass_variants'] += 1

            # Get VEP annotations
            csq = variant.INFO.get('CSQ')
            if not csq:
                continue

            consequences = self.parse_vep_consequences(csq)

            # Count consequences for all variants (split by coding/non-coding)
            for cons in consequences:
                if cons in CODING_CONSEQUENCES:
                    self.stats['coding_consequence_counts'][cons] = self.stats['coding_consequence_counts'].get(cons, 0) + 1
                else:
                    self.stats['noncoding_consequence_counts'][cons] = self.stats['noncoding_consequence_counts'].get(cons, 0) + 1

            # Classify variant types
            is_coding_variant = self.is_coding(consequences)
            is_nonsynonymous_variant = self.is_nonsynonymous(consequences)
            is_synonymous_variant = self.is_synonymous(consequences)

            if is_coding_variant:
                self.stats['coding_variants'] += 1
            else:
                self.stats['noncoding_variants'] += 1

            if is_nonsynonymous_variant:
                self.stats['nonsynonymous_variants'] += 1

            if is_synonymous_variant:
                self.stats['synonymous_variants'] += 1

            # Apply depth and VAF filters to ALL variants
            dp_arr = variant.format('DP')
            vaf_arr = variant.format('VAF')

            depth_ok = True
            vaf_ok = True
            depth = None
            variant_af = None

            # Both of these assume there's only a single value per variant (single-sample enforced)
            if dp_arr is not None and len(dp_arr) > 0:
                depth = dp_arr[0]
                try:
                    if hasattr(depth, '__len__') and not isinstance(depth, (str, bytes)):
                        depth = depth[0]
                except Exception:
                    pass
                if depth is not None and depth < self.min_depth:
                    depth_ok = False

            if vaf_arr is not None and len(vaf_arr) > 0:
                variant_af = vaf_arr[0]
                try:
                    if hasattr(variant_af, '__len__') and not isinstance(variant_af, (str, bytes)):
                        variant_af = variant_af[0]
                except Exception:
                    pass
                if variant_af is not None and variant_af < self.min_vaf:
                    vaf_ok = False

            # Compute gnomAD MAX_AF from VEP CSQ (null/empty -> 0.0)
            max_af = self.vep_parser.get_max_af(csq)
            pass_gnomAD = max_af <= self.max_af_threshold

            # Record debug info before early-continue if requested
            if self.debug_tsv:
                self._debug_rows.append({
                    'CHROM': variant.CHROM,
                    'POS': variant.POS,
                    'REF': variant.REF,
                    'ALT': ','.join(variant.ALT or []),
                    'Consequence': '&'.join(sorted(set(consequences))) if consequences else '',
                    'MAX_AF': max_af,
                    'DP': depth if depth is not None else '',
                    'VAF': float(variant_af) if variant_af is not None else '',
                    'IsCoding': is_coding_variant,
                    'IsNonsynonymous': is_nonsynonymous_variant,
                    'PassDepthVAF': bool(depth_ok and vaf_ok),
                    'PassGnomAD': bool(pass_gnomAD),
                    'CountedInTMB': False
                })

            if not (depth_ok and vaf_ok):
                continue

            if not pass_gnomAD:
                # Exclude common variants
                self.stats.setdefault('gnomad_filtered_variants', 0)
                self.stats['gnomad_filtered_variants'] += 1
                continue

            # Count consequences for filtered variants (split by coding/non-coding)
            for cons in consequences:
                if cons in CODING_CONSEQUENCES:
                    self.stats['filtered_coding_consequence_counts'][cons] = self.stats['filtered_coding_consequence_counts'].get(cons, 0) + 1
                else:
                    self.stats['filtered_noncoding_consequence_counts'][cons] = self.stats['filtered_noncoding_consequence_counts'].get(cons, 0) + 1

            # Count ALL qualifying variants for TMB
            qualifying_all_variants += 1
            self.stats['filtered_all_variants'] += 1
            if self.debug_tsv:
                self._debug_rows[-1]['CountedInTMB'] = True

            # Count qualifying coding variants separately
            if is_coding_variant:
                qualifying_coding_variants += 1
                self.stats['filtered_coding_variants'] += 1

            # Count qualifying nonsynonymous variants separately
            if is_nonsynonymous_variant:
                qualifying_nonsynonymous_variants += 1
                self.stats['filtered_nonsynonymous_variants'] += 1

        console.print(f"[green]Found {qualifying_all_variants:,} total qualifying variants (for TMB)")
        console.print(f"[green]Found {qualifying_coding_variants:,} qualifying coding variants")
        console.print(f"[green]Found {qualifying_nonsynonymous_variants:,} qualifying nonsynonymous variants")
        return qualifying_all_variants, qualifying_coding_variants, qualifying_nonsynonymous_variants

    def calculate_tmb(self) -> Dict:
        """Calculate TMB score and generate report"""
        console.print("[bold blue]Calculating TMB...")

        # Calculate eligible regions
        eligible_bases, avg_coverage = self.calculate_eligible_regions()
        eligible_mbp = eligible_bases / 1_000_000

        # Process variants
        qualifying_all_variants, qualifying_coding_variants, qualifying_nonsynonymous_variants = self.process_variants()

        # Calculate TMB (based on ALL variants - coding and non-coding)
        if eligible_mbp > 0:
            tmb_score = qualifying_all_variants / eligible_mbp
        else:
            tmb_score = 0
            console.print("[red]Warning: No eligible regions found!")

        # Denominator for nonsynonymous score
        # If no region BED provided, assume CDS size is 34 Mbp and VCF is whole genome
        if not self._bed_intervals:
            nonsyn_den_mbp = 34.0
            self.stats['assumed_cds_mbp'] = 34.0
            self.stats['assumed_vcf_whole_genome'] = True
        else:
            nonsyn_den_mbp = eligible_mbp

        # Calculate nonsynonymous score using selected denominator
        if nonsyn_den_mbp and nonsyn_den_mbp > 0:
            nonsynonymous_score = qualifying_nonsynonymous_variants / nonsyn_den_mbp
        else:
            nonsynonymous_score = 0

        # Record denominator used for transparency
        self.stats['nonsynonymous_denominator_mbp'] = float(nonsyn_den_mbp)

        self.stats['tmb_score'] = tmb_score
        self.stats['nonsynonymous_score'] = nonsynonymous_score

        # Optionally write debug TSV
        if self.debug_tsv and self._debug_rows:
            self._write_debug_tsv(self.debug_tsv)

        # Display results
        self.display_results(show_consequences=self.show_consequences)

        return self.stats

    def display_consequence_table(self, show_consequences: bool = False):
        """Display consequence counts tables if requested"""
        if not show_consequences:
            return

        # Display coding consequences
        if self.stats['coding_consequence_counts']:
            console.print()
            table = Table(title="Coding Variant Consequences", show_header=True, header_style="bold blue")
            table.add_column("Consequence", style="cyan")
            table.add_column("All Variants", style="white", justify="right")
            table.add_column("Filtered Variants", style="green", justify="right")

            # Sort consequences by frequency (all variants)
            sorted_consequences = sorted(self.stats['coding_consequence_counts'].items(),
                                       key=lambda x: x[1], reverse=True)

            for consequence, all_count in sorted_consequences:
                filtered_count = self.stats['filtered_coding_consequence_counts'].get(consequence, 0)
                table.add_row(
                    consequence.replace('_', ' ').title(),
                    f"{all_count:,}",
                    f"{filtered_count:,}"
                )

            console.print(table)

        # Display non-coding consequences
        if self.stats['noncoding_consequence_counts']:
            console.print()
            table = Table(title="Non-coding Variant Consequences", show_header=True, header_style="bold yellow")
            table.add_column("Consequence", style="cyan")
            table.add_column("All Variants", style="white", justify="right")
            table.add_column("Filtered Variants", style="green", justify="right")

            # Sort consequences by frequency (all variants)
            sorted_consequences = sorted(self.stats['noncoding_consequence_counts'].items(),
                                       key=lambda x: x[1], reverse=True)

            for consequence, all_count in sorted_consequences:
                filtered_count = self.stats['filtered_noncoding_consequence_counts'].get(consequence, 0)
                table.add_row(
                    consequence.replace('_', ' ').title(),
                    f"{all_count:,}",
                    f"{filtered_count:,}"
                )

            console.print(table)

    def display_results(self, show_consequences: bool = False):
        """Display TMB results in a formatted table"""
        table = Table(title="TMB Calculation Results", show_header=True, header_style="bold magenta")
        table.add_column("Metric", style="cyan", no_wrap=True)
        table.add_column("Value", style="magenta")

        table.add_row("Total variants", f"{self.stats['total_variants']:,}")
        table.add_row("PASS variants", f"{self.stats['pass_variants']:,}")
        table.add_row("├─ Coding variants", f"{self.stats['coding_variants']:,}")
        table.add_row("│  ├─ Nonsynonymous", f"{self.stats['nonsynonymous_variants']:,}")
        table.add_row("│  └─ Synonymous", f"{self.stats['synonymous_variants']:,}")
        table.add_row("└─ Non-coding variants", f"{self.stats['noncoding_variants']:,}")
        table.add_row("", "")
        table.add_row("Qualifying ALL variants", f"{self.stats['filtered_all_variants']:,}")
        table.add_row("├─ Qualifying coding", f"{self.stats['filtered_coding_variants']:,}")
        table.add_row("└─ Qualifying nonsynonymous", f"{self.stats['filtered_nonsynonymous_variants']:,}")
        table.add_row("", "")
        table.add_row("Average coverage", f"{self.stats['average_coverage']:.2f}x")
        table.add_row("Eligible bases", f"{self.stats['eligible_bases']:,}")
        table.add_row("Eligible region", f"{self.stats['eligible_mbp']:.2f} Mbp")
        # Explicit note when region BED used
        if self.stats.get('region_bed_used'):
            table.add_row("Region BED used", "Yes")
            # Show size of the BED for context
            table.add_row("Region BED size", f"{self.stats.get('region_bed_total_mbp', 0.0):.2f} Mbp")
        else:
            # When no region BED is provided, we assume whole-genome VCF and CDS ~34 Mbp
            table.add_row("Assumption: VCF scope", "Whole genome")
            table.add_row("Assumption: CDS size", f"{self.stats.get('assumed_cds_mbp', 34.0):.2f} Mbp")
            table.add_row("Nonsyn denom (Mbp)", f"{self.stats.get('nonsynonymous_denominator_mbp', 34.0):.2f}")
        table.add_row("", "")
        table.add_row("[bold]TMB Score (ALL Variants)[/bold]", f"[bold green]{self.stats['tmb_score']:.2f} mutations/Mbp[/bold green]")
        table.add_row("[bold]Nonsynonymous Score[/bold]", f"[bold yellow]{self.stats['nonsynonymous_score']:.2f} mutations/Mbp[/bold yellow]")
        if 'gnomad_filtered_variants' in self.stats:
            table.add_row("", "")
            table.add_row(
                f"gnomAD MAX_AF > {self.max_af_threshold*100:.2f}% filtered",
                f"{self.stats['gnomad_filtered_variants']:,}"
            )

        console.print(table)

        # Display consequence breakdown if requested
        self.display_consequence_table(show_consequences)

    def _write_debug_tsv(self, path: str):
        """Write collected per-variant debug information to a TSV file."""
        fieldnames = [
            'CHROM', 'POS', 'REF', 'ALT',
            'Consequence', 'MAX_AF', 'DP', 'VAF',
            'IsCoding', 'IsNonsynonymous',
            'PassDepthVAF', 'PassGnomAD', 'CountedInTMB'
        ]
        try:
            with open(path, 'w', newline='') as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                for row in self._debug_rows:
                    writer.writerow(row)
            console.print(f"[green]Debug TSV written to {path}")
        except Exception as e:
            console.print(f"[red]Failed to write debug TSV to {path}: {e}")


@click.command()
@click.option('--vcf', required=True, help='VEP-annotated VCF file')
@click.option('--coverage', required=True, help='Mosdepth coverage BED file')
@click.option(
    '--region-bed',
    type=str,
    help=(
        'BED file of regions to restrict TMB (e.g., coding regions). '
        'If omitted, the nonsynonymous score uses a constant 34 Mbp denominator '
        '(assumes whole-genome VCF and CDS size ~34 Mbp).'
    ),
)
@click.option('--min-vaf', default=0.10, type=float, help='Minimum variant allele frequency (default: 0.10)')
@click.option('--min-depth', default=10, type=int, help='Minimum read depth (default: 10)')
@click.option('--min-coverage', default=10, type=int, help='Minimum coverage for eligible regions (default: 10)')
@click.option('--max-af-cutoff', default=0.03, type=float, show_default=True,
              help='Exclude variants with gnomAD MAX_AF greater than this value')
@click.option('--show-consequences', is_flag=True, help='Display detailed consequence breakdown table')
@click.option('--debug-tsv', type=str, help='Write per-variant debug TSV with VEP fields used')
@click.option('--output', help='Output JSON file for results')
def main(vcf, coverage, region_bed, min_vaf, min_depth, min_coverage, max_af_cutoff, show_consequences, debug_tsv, output):
    """Calculate Tumor Mutational Burden (TMB) from VEP-annotated VCF and coverage files.

    Notes:
    - TMB score (ALL variants) uses eligible region derived from coverage and optional --region-bed.
    - Nonsynonymous score denominator:
      * With --region-bed: uses eligible region (Mbp).
      * Without --region-bed: uses a constant 34 Mbp (assumes whole-genome VCF and CDS size ~34 Mbp).
    """

    # Check input files exist
    if not Path(vcf).exists():
        console.print(f"[red]Error: VCF file not found: {vcf}")
        sys.exit(1)

    if not Path(coverage).exists():
        console.print(f"[red]Error: Coverage file not found: {coverage}")
        sys.exit(1)

    console.print(f"[bold green]TMB Calculator[/bold green]")
    console.print(f"VCF file: {vcf}")
    console.print(f"Coverage file: {coverage}")
    if region_bed:
        console.print(f"Region BED: {region_bed}")
    console.print(f"Filters: VAF ≥ {min_vaf}, Depth ≥ {min_depth}, Coverage ≥ {min_coverage}, Exclude MAX_AF > {max_af_cutoff}")
    console.print()

    # Initialize calculator
    calculator = TMBCalculator(
        vcf_path=vcf,
        coverage_path=coverage,
        min_vaf=min_vaf,
        min_depth=min_depth,
        min_coverage=min_coverage,
        show_consequences=show_consequences,
        max_af_threshold=max_af_cutoff,
        debug_tsv=debug_tsv,
        region_bed=region_bed
    )

    # Calculate TMB
    results = calculator.calculate_tmb()

    # Save results if output specified
    if output:
        with open(output, 'w') as f:
            json.dump(results, f, indent=2)
        console.print(f"[green]Results saved to {output}")


if __name__ == '__main__':
    main()
