# TMB Calculator

Tumor Mutational Burden (TMB) calculator for VEP-annotated VCF files

## Overview

This tool calculates TMB from somatic variants in VEP-annotated VCF files, using coverage information from mosdepth to determine eligible regions. TMB is calculated using **all variants** (coding and non-coding) as the numerator, with a separate score for nonsynonymous variants. Variants can be filtered by gnomAD frequency via the VEP `MAX_AF` annotation. Optionally, limit both the eligible region and counted variants to a provided BED (e.g., coding regions inferred from GENCODE GTF). If no region BED is provided, the nonsynonymous score uses a constant 34 Mbp denominator (assuming whole-genome VCF and CDS size of 34 Mbp) and the console/JSON output explicitly states this assumption.

## Features

- **Correct TMB calculation**: Uses ALL variants (coding + non-coding) as numerator
- **Separate nonsynonymous score**: Additional metric for protein-altering variants only
- **34 Mbp fallback**: Without a region BED, nonsynonymous score uses a constant 34 Mbp denominator and prints the assumption in console and JSON
- **VEP annotation parsing**: Extracts fields from the CSQ header format
- **gnomAD filtering**: Excludes variants with VEP `MAX_AF` above a configurable cutoff
- **Consequence breakdown**: Optional detailed tables showing coding and non-coding variant counts by consequence type
- **Configurable filtering**: Adjustable VAF, depth, and coverage thresholds
- **Rich console output**: Progress tracking and formatted results
- **JSON export**: Machine-readable results for downstream analysis

## Installation

Using pixi (recommended):

```bash
pixi install
```

Example with gnomAD filtering, coding regions BED, and debug TSV:

```bash
pixi run python calculate_tmb.py \
  --vcf your_variants.vep.vcf.gz \
  --coverage your_coverage.regions.bed.gz \
  --region-bed coding_regions.bed.gz \
  --min-vaf 0.10 \
  --min-depth 10 \
  --min-coverage 10 \
  --max-af-cutoff 0.03 \
  --debug-tsv variants_debug.tsv \
  --output results.json
```

## Usage

### Basic Usage

```bash
# Test with COLO829 data
pixi run test

# Test with detailed consequence breakdown
pixi run test-detailed

# Custom parameters
pixi run python calculate_tmb.py \
  --vcf your_variants.vep.vcf.gz \
  --coverage your_coverage.regions.bed.gz \
  --min-vaf 0.10 \
  --min-depth 10 \
  --min-coverage 10 \
  --max-af-cutoff 0.03 \
  --debug-tsv variants_debug.tsv \
  --show-consequences \
  --output results.json
```

### Available Scripts

1. **`calculate_tmb.py`** (Main script): TMB calculation using all variants with VEP parsing and optional gnomAD `MAX_AF` filtering; optionally restrict to a BED with `--region-bed`
2. **`gtf_to_coding_bed.py`**: Convert a GENCODE GTF (e.g., `gencode.v46.annotation.gtf.gz`) to a merged BED of coding regions (protein-coding CDS only), outputting a 4-column BED with gene name

### Input Files

- **VCF file**: VEP-annotated somatic variants (must contain CSQ field)
  - Must be a single-sample VCF; the tool fails if more than one sample is present
  - Must include a VEP CSQ INFO header with a `Format:` list containing `Consequence` and `MAX_AF` (required)
  - If `MAX_AF` is missing/empty (`.`), it is treated as 0.0
  - Assumes single values per variant for depth (`DP`) and allele frequency (`VAF`) fields
- **Coverage file**: Mosdepth regions.bed.gz with 4 columns (chr, start, end, coverage)

### Parameters

- `--min-vaf`: Minimum variant allele frequency (default: 0.10)
- `--min-depth`: Minimum read depth (default: 10) 
- `--min-coverage`: Minimum coverage for eligible regions (default: 10)
- `--show-consequences`: Display detailed consequence breakdown table (flag)
- `--max-af-cutoff`: Exclude variants with VEP gnomAD `MAX_AF` greater than this value (default: 0.03)
- `--debug-tsv`: Write a per-variant TSV with relevant fields used during filtering
- `--region-bed`: BED file to restrict both eligible region (denominator) and counted variants (numerator); e.g., coding regions BED from `gtf_to_coding_bed.py`
  - BED may have extra columns; only the first three are used for intersection.

### Generating a Coding Regions BED

Convert a GENCODE GTF to a BED of coding regions (merged CDS intervals for protein-coding transcripts, per gene; includes gene name in 4th column):

```bash
python gtf_to_coding_bed.py gencode.v46.annotation.gtf.gz -o coding_regions.bed.gz
```

Notes:
- Uses `feature == CDS` and `transcript_type == protein_coding` (or `transcript_biotype`).
- Includes primary chromosomes `chr1..chr22`, `chrX`, `chrY` by default; pass `--include-chrM` to include chrM.
- Output is 0-based half-open BED, merged per gene, with the 4th column set to `gene_name`.

### Output Metadata When Using a BED

When `--region-bed` is provided, the console table and the JSON include BED metadata:
- `region_bed_used`: whether a BED was used
- `region_bed_path`: path to the BED
- `region_bed_total_bases` and `region_bed_total_mbp`: size of the BED after chrM filtering

## TMB Calculation Details

### TMB Formula

```
TMB (ALL variants) = ALL Qualifying Variants / Eligible Region (Mbp)

Nonsynonymous Score =
  - With --region-bed: Nonsynonymous Variants / Eligible Region (Mbp)
  - Without --region-bed: Nonsynonymous Variants / 34 Mbp
```

Notes:
- When no region BED is provided, we assume the VCF represents whole-genome calls and the total CDS size is 34 Mbp. The tool prints these assumptions in the console and includes them in the JSON.

### Variant Classification

**Nonsynonymous consequences** (for separate score):
- missense_variant
- frameshift_variant
- stop_gained, stop_lost, start_lost
- splice_acceptor_variant, splice_donor_variant
- inframe_insertion, inframe_deletion
- protein_altering_variant
- feature_elongation, feature_truncation
- incomplete_terminal_codon_variant, transcript_truncation

**Coding consequences** (includes synonymous):
- All nonsynonymous consequences above
- synonymous_variant
- stop_retained_variant, start_retained_variant

### Filtering Criteria

1. **Variant quality**: Only PASS variants
2. **Genomic regions**: Excludes mitochondrial variants
3. **Coverage**: Only regions with coverage ≥ threshold
4. **Variant properties**: VAF and depth thresholds (single-sample format assumed)
5. **Population frequency**: Exclude variants with VEP `MAX_AF` > `--max-af-cutoff` (requires `MAX_AF` in CSQ header; missing values treated as 0.0)

### Format Assumptions

- **VCF Format**: DeepSomatic output with VEP annotations
- **Single Sample**: Exactly one sample is required; the tool fails if multiple samples are present
- **Field Access**: Assumes single values for DP and VAF format fields per variant
- **VEP CSQ**: Must contain a `##INFO=<ID=CSQ,...>` header with a `Format:` list including `Consequence` and `MAX_AF`
- **VEP Annotations**: Based on one consequence per transcript (using `--pick` flag)
- **VEP Parameters**: For full VEP annotation parameters, refer to the hifi-somatic-wdl pipeline

## Test Results (COLO829)

```
TMB Calculation Results
├── Total variants: 42,662
├── ├─ Coding variants: 351
├── │  ├─ Nonsynonymous: 225
├── │  └─ Synonymous: 126
├── └─ Non-coding variants: 42,311
├── Qualifying ALL variants: 42,316
├── ├─ Qualifying coding: 345
├── └─ Qualifying nonsynonymous: 221
├── Eligible region: 2,882.21 Mbp
├── TMB Score (ALL Variants): 14.68 mutations/Mbp
└── Nonsynonymous Score: 0.08 mutations/Mbp
```

### Example Output When No Region BED Is Provided

The console table includes explicit assumption lines and the denominator used for the nonsynonymous score:

```
TMB Calculation Results
...
Eligible region: 2,882.21 Mbp
Assumption: VCF scope: Whole genome
Assumption: CDS size: 34.00 Mbp
Nonsyn denom (Mbp): 34.00

TMB Score (ALL Variants): 14.68 mutations/Mbp
Nonsynonymous Score: 6.50 mutations/Mbp
```

**Key Insights:**
- TMB (all variants): 14.68 mut/Mbp - close to expected ~10 mut/Mbp for COLO829
- Nonsynonymous score: 0.08 mut/Mbp - much lower, as expected
- Most variants are non-coding (intron/intergenic), with coding variants representing <1% of total

## Interpretation

TMB interpretation typically uses the **all variants** score:
- **High TMB** (≥10 mut/Mbp): May benefit from immunotherapy
- **Intermediate TMB** (5-10 mut/Mbp): Variable clinical significance  
- **Low TMB** (<5 mut/Mbp): Less likely to respond to immunotherapy

## Output Format

Results are saved as JSON:

```json
{
  "total_variants": 42662,
  "coding_variants": 351,
  "nonsynonymous_variants": 225,
  "synonymous_variants": 126,
  "noncoding_variants": 42311,
  "filtered_all_variants": 42316,
  "filtered_coding_variants": 345,
  "filtered_nonsynonymous_variants": 221,
  "eligible_mbp": 2882.21,
  "tmb_score": 14.68,
  "nonsynonymous_score": 0.08,
  "average_coverage": 63.04,
  "coding_consequence_counts": {
    "missense_variant": 194,
    "synonymous_variant": 126,
    "stop_gained": 13,
    "frameshift_variant": 7
  },
  "noncoding_consequence_counts": {
    "intron_variant": 18866,
    "intergenic_variant": 17873,
    "upstream_gene_variant": 1675,
    "regulatory_region_variant": 1617
  },
  "filtered_coding_consequence_counts": {
    "missense_variant": 190,
    "synonymous_variant": 124,
    "stop_gained": 13,
    "frameshift_variant": 7
  },
  "filtered_noncoding_consequence_counts": {
    "intron_variant": 18707,
    "intergenic_variant": 17737,
    "upstream_gene_variant": 1654,
    "regulatory_region_variant": 1603
  }
}
```

When `--max-af-cutoff` is set, the console summary includes a line showing how many variants were filtered due to exceeding the cutoff. If `--debug-tsv` is provided, a TSV file is written with columns: `CHROM, POS, REF, ALT, Consequence, MAX_AF, DP, VAF, IsCoding, IsNonsynonymous, PassDepthVAF, PassGnomAD, CountedInTMB`.

When no region BED is provided, additional JSON fields capture the assumptions and denominator used for the nonsynonymous score:

```json
{
  "assumed_cds_mbp": 34.0,
  "assumed_vcf_whole_genome": true,
  "nonsynonymous_denominator_mbp": 34.0
}
```

Field meanings:
- `assumed_cds_mbp`: Set to 34.0 when `--region-bed` is not used (otherwise null).
- `assumed_vcf_whole_genome`: true when no `--region-bed` is used; indicates WGS scope assumption.
- `nonsynonymous_denominator_mbp`: The denominator used for the nonsynonymous score: `eligible_mbp` with a BED, or 34.0 without a BED.

## Dependencies

- cyvcf2: VCF parsing
- rich: Console output formatting  
- click: Command-line interface
- numpy, pandas: Data processing
