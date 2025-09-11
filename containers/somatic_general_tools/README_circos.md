# Circos Plot Generator for Structural Variants

A Python script to generate circos plots from VCF files containing structural variants, with intelligent label positioning and overlap prevention.

## Requirements

- Python environment with pycirclize
- Required files: `hg38.bed`, `hg38_cytoband.tsv`, Mitelman database file

## Basic Usage

```bash
python plot_circos_2.py <vcf_file> <sample_name> <mitelman_mcgene_file> [options]
```

### Required Arguments

- `vcf_file`: Input VCF file (.vcf or .vcf.gz)
- `sample_name`: Sample name for output files
- `mitelman_mcgene_file`: Mitelman fusion database file (e.g., MCGENE.TXT.DATA)

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--max-gene-labels` | 60 | Maximum gene labels to display (0 = no limit) |
| `--chr-label-size` | 10 | Chromosome label font size |
| `--sector-space` | 5 | Space between chromosomes |
| `--link-alpha` | 0.8 | Translocation line transparency (0.0-1.0) |

## Examples

### Default (Recommended)
```bash
python plot_circos_2.py sample.vcf.gz SAMPLE_001 MCGENE.TXT.DATA
```

### Clean Plot for Complex Data
```bash
python plot_circos_2.py sample.vcf.gz SAMPLE_001 MCGENE.TXT.DATA \
  --sector-space 8 \
  --link-alpha 0.6 \
  --max-gene-labels 40
```

### Increased Spacing for Better Readability
```bash
python plot_circos_2.py sample.vcf.gz SAMPLE_001 MCGENE.TXT.DATA \
  --sector-space 10
```

### Show All Genes (No Limit)
```bash
python plot_circos_2.py sample.vcf.gz SAMPLE_001 MCGENE.TXT.DATA \
  --max-gene-labels 0
```

## Output Files

- `<sample_name>.pdf`: High-resolution PDF
- `<sample_name>.png`: High-resolution PNG

## Features

- **Clean Chromosome Labels**: Labels positioned below chromosome segments without overlapping
- **Known Fusion Highlighting**: Red lines and text for known fusion pairs
- **Priority System**: Known fusions labeled first when gene limit is reached
- **Smart Gene Grouping**: Nearby genes combined to prevent overlap
- **Transparency**: Reduces visual clutter from overlapping lines
- **Clean Naming**: Removes "chr" prefix from chromosome labels (1, 2, X, Y instead of chr1, chr2, etc.)

## Tips

- Increase `--sector-space` (8-12) for better chromosome separation and readability
- Reduce `--link-alpha` (0.5-0.7) to make translocation lines less prominent
- Set `--max-gene-labels 0` only for datasets with few structural variants
- For complex datasets with many translocations, use increased spacing (10+)