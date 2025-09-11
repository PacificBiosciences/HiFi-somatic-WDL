# Output from the workflow

- [Output from the workflow](#output-from-the-workflow)
  - [QC results](#qc-results)
  - [Structural variants](#structural-variants)
  - [SNV/INDEL (Somatic and germline)](#snvindel-somatic-and-germline)
  - [Phasing](#phasing)
  - [Homologous-recombination deficiency prediction](#homologous-recombination-deficiency-prediction)
  - [Differentially methylated region](#differentially-methylated-region)
  - [Somatic SNV/INDEL annotation](#somatic-snvindel-annotation)
  - [Prioritization](#prioritization)
  - [Mutational signatures](#mutational-signatures)
  - [Purity and ploidy estimation](#purity-and-ploidy-estimation)
  - [Copy number variation](#copy-number-variation)
  - [Biomarkers](#biomarkers)
  - [Differentially methylated regions](#differentially-methylated-regions)
  - [Report HTML](#report-html)

## QC results

Depth of coverage for both tumor and normal samples is in the `mosdepth_normal_summary` and `mosdepth_tumor_summary` folders. The `mosdepth` folder contains per-chromosome depth and overall coverage (row `total`).

`overall_(tumor\|normal)_alignment_stats` contains alignment statistics for both tumor and normal samples, generated using `seqkit bam`. Per-alignment statistics are in `per_alignment_(tumor\|normal)_stats`.

## Structural variants

The pipeline incorporates the structural variant caller [Severus](https://github.com/KolmogorovLab/Severus). Results are filtered using a germline SV VCF (Human Pangenome Reference Consortium) via `svpack` to reduce false positives. `svpack` also adds basic annotations from the Ensembl GFF (v101). The filtered VCFs are in `Severus_filtered_vcf`.

Severus can visualize complex SV clusters. All cluster HTMLs are archived as `Severus_cluster_plots`.

`AnnotSV` is used to further annotate SVs to a TSV file (`AnnotatedSeverusSV`). The TSV format is described in the [AnnotSV README](https://github.com/lgmgeo/AnnotSV/blob/master/README.AnnotSV_latest.pdf). This provides more detailed annotation than `svpack`.

Finally, to help prioritize cancer-relevant variants, the workflow annotates SVs with the IntOGen Compendium of Cancer Genes (CCG) and produces the final set in `Annotated*SV_intogen`.

In `SV_circos`, we generate a circos plot of all BND events at least 100 kbp apart, useful for visualizing large-scale rearrangements. Red denotes fusion pairs recorded in the Mitelman database that has more than 3 records (samples); black denotes fusions involving genes in the `svpack` GFF; grey denotes fusions not in either source. We also output a list of breakpoints pairs at `SV_known_genes_pairs` when at least one of the genes is found in Mitelman database (3 samples).

## SNV/INDEL (Somatic and germline)

Somatic SNV/INDEL VCFs are in `small_variant_vcf`. We use DeepSomatic (default) or ClairS to call somatic SNV/INDELs. See [input JSON parameters](step-by-step.md#input-json-parameters) to switch callers or select the Sequel IIe ClairS model.

DeepSomatic is currently computationally expensive (14–18 hours for 60X/30X tumor/normal) and requires a separate caller for germline SNV/INDELs. The workflow uses Clair3 to call germline variants in both tumor and normal. DeepSomatic does not report normal depth at variant sites, so its VCF cannot be used as input to Purple for purity and ploidy estimation. If DeepSomatic is used, Purple runs without a somatic VCF, which may affect estimates.

For both ClairS and DeepSomatic, the workflow can split the genome into chunks (default 75 Mbp) and call variants in parallel, then merge outputs into a single VCF. This scales to large genomes and datasets by leveraging multiple HPC nodes. Germline variants are called with Clair3 regardless of the somatic caller.

For annotation, we use Ensembl VEP to annotate the VCFs (see `small_variant_tsv_annotated` and `small_variant_vcf_annotated`). As with SVs, the workflow also annotates SNV/INDELs with the IntOGen CCG and produces the final set in `small_variant_tsv_CCG`.

## Phasing

The workflow uses HiPhase to phase both germline and somatic variants. As a diploid phasing tool, HiPhase forces somatic variants into parental haplotypes rather than reflecting tumor polyclonality, but this is still useful for understanding haplotype structure. In the final VCF, the `PS` tag represents phase blocks (i.e., variants with the same `PS` tag belong to the same haplotype within a block). We set all somatic genotypes in DeepSomatic output to `0/1` (DeepSomatic defaults to `1/1`) to enable HiPhase phasing.

## Homologous-recombination deficiency prediction

SNVs, INDELs, and SVs are supplied to CHORD for HRD prediction. Results are in `chord_hrd_prediction`. Internal testing on two HRD-positive datasets produced accurate results. However, broader validation is pending, so use with caution—especially at low tumor purity and coverage. For example, at 15X effective tumor coverage (30X at 50% purity or 60X at 25%), CHORD predicted HCC1395 to be HR-deficient, which is not true.

## Differentially methylated regions

CpG calls at each locus in the human genome are summarized using [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools). The BED file of CpG calls is then used to call DMRs with DSS. The pileup BED file is in `pileup_(normal\/tumor)_bed`.

`DMR_annotated` contains five folders per patient. Each folder contains DMRs annotated to different genomic features. Annotation is done using [annotatr](https://bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html).

``` markdown
├── 0
│   └── COLO829_hg38_genes_1to5kb_dmrs_intogen.nCG50_summary.tsv.gz -> DMRs that are 1 to 5 kbp away from TSS of genes
├── 1
│   └── COLO829_hg38_genes_3UTRs_dmrs_intogen.nCG50_summary.tsv.gz ->  DMRs that are 3'UTR of genes
├── 2
│   └── COLO829_hg38_genes_5UTRs_dmrs_intogen.nCG50_summary.tsv.gz ->  DMRs that are 5'UTR of genes
├── 3
│   └── COLO829_hg38_genes_exons_dmrs_intogen.nCG50_summary.tsv.gz -> DMRs that are exons of genes
├── 4
│   └── COLO829_hg38_genes_introns_dmrs_intogen.nCG50_summary.tsv.gz -> DMRs that are introns of genes
└── 5
    └── COLO829_hg38_genes_promoters_dmrs_intogen.nCG50_summary.tsv.gz -> DMRs that are (known) promoters of genes
```

For each TSV, the first 11 columns (ending with `areaStat`) are identical and are produced by [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html).
- `meanMethyl1`: mean methylation level in tumor.
- `meanMethyl2`: mean methylation level in normal.
- `length`: length of the DMR.
- `nCG`: number of CpG sites in the DMR. By default, the workflow requires at least 50 CpG sites per DMR.
- `areaStat`: DMR area statistic; larger values indicate greater significance.

`annot.*` columns are produced by `annotatr`; all uppercase columns are extracted from the IntOGen Compendium of Cancer Genes TSV file.

## Somatic SNV/INDEL annotation

`small_variant_vcf_annotated` contains somatic variants annotated using Ensembl VEP. To obtain readable annotations, try the `split-vep` tool from [`bcftools`](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html). The workflow also produces a TSV using `split-vep` in `small_variant_tsv`.

## Prioritization

The workflow prioritizes variants using the IntOGen CCG: DMRs overlapping CCG genes with ≥50 CpGs, and SV/SNV/INDELs involving any CCG gene.

## Mutational signatures

Mutational signatures are determined using the `MutationalPatterns` R [package](https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html).
There are two folders in the output:
- `mutsig_SNV` contains TSV files of mutational signatures fitted to known COSMIC signatures. The TSVs are:
  - `mut_sigs.tsv`: contribution of the fitted signatures (using `fit_to_signatures_strict`).
  - `reconstructed_sigs.tsv`: reconstructed trinucleotide context based on the fitted signatures.
  - `type_occurences.tsv`: frequency of each mutation type in the sample.
  - `mut_sigs_bootstrapped.tsv`: bootstrapped fit results (via `fit_to_signatures_bootstrapped`), indicating signature stability.

## Purity and ploidy estimation

The workflow implements the HMFtools suite to estimate purity and ploidy using Amber, Cobalt, and Purple (tumor/normal workflow). 
For Cobalt, due to noise in long-read-based read depth and B-allele frequency, we set PCF gamma to 1000 for better segmentation. See the [GitHub issue](https://github.com/hartwigmedical/hmftools/issues/485) for discussion. In practice, CNV calls may miss smaller focal events, but those should be picked up by the SV caller.

Purity and ploidy estimates have been reasonably robust in our experience with HCC1395 and COLO829, 
but should be used with caution. Files ending with `.purity.tsv` are in the `Purple_outputs` folder.

Wakhan also estimates purity and ploidy; see the section below on copy number variation.

## Copy number variation

CNVKit segments copy numbers from matched tumor/normal BAMs. For long reads, we use 10 kbp bins, which performed best on COLO829. The workflow also uses HMFtools purity/ploidy estimates together with Clair's heterozygous SNVs to estimate allele-specific major and minor copy numbers in `cnvkit_cns_with_major_minor_CN`.

CNVKit's recall mode requires an integer ploidy, which can fail with subclonal CNV.
Purple also calls allele-specific copy numbers and can account for non-integer ploidy (subclonal CNV). Purple CNV segments are in `Purple_outputs` with the suffix `tumor.purple.cnv.somatic.tsv`. We visualize Purple results in the final report.

Lastly, we include [Wakhan](https://github.com/KolmogorovLab/Wakhan) as an alternative CNV caller for both tumor/normal and tumor-only workflows. Wakhan calls haplotype-specific CNV and also estimates purity and ploidy. Results are in `wakhan_cnv`. We select the best purity and ploidy solution in the final results folder; the full set of solutions is available by extracting the `tar.gz` file. As of v0.9, Wakhan is still in active development, so we report Purple by default for tumor/normal. For tumor-only, we use Wakhan (Purple has not been tested for tumor-only).

## Biomarkers

We developed `owl` to profile microsatellite repeats in tumors (no normal required). Each site is scanned for per-read length, and the variability at each site is summarized in `owl_msi_profile`. A score is calculated based on the percentage of unstable sites (`owl_msi_score`). Based on internal data across MSI-high tumors and background populations, >10% unstable sites is a good indicator of MSI.

We also provide a Python script to estimate tumor mutation burden (TMB). It takes the annotated small variants as input and categorizes mutations into synonymous and non-synonymous. If a region BED file is supplied, only mutations within those regions are counted. By default, the pipeline reports a whole-genome TMB and a non-synonymous TMB based on an estimated 35 Mbp of coding sequence in the human genome. It also reports an estimate using Gencode v46 coding regions.

## Report HTML

The workflow produces an HTML report summarizing the results, available in the `report` folder. It also includes variants filtered by documented criteria to aid interpretation. These criteria are guidance, not strict rules, and may not be optimal for all datasets. We encourage users to analyze the pipeline outputs in more detail.
