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
  - [Report html](#report-html)

## QC results

Depth of coverage for both tumor and normal can be found in the `mosdepth_normal_summary` and `mosdepth_tumor_summary` folders. The `mosdepth` folder contains the depth of coverage for each chromosome and overall coverage (`total` at the end of the file).

`overall_(tumor\|normal)_alignment_stats` contains the alignment statistics for both tumor and normal samples. The statistics are produced by using the `seqkit bam` command. Per-alignment statistics are also generated in the `per_alignment_(tumor\|normal)_stats` folders.

## Structural variants

The pipeline incorporates the structural variants caller [Severus](https://github.com/KolmogorovLab/Severus) with default parameters. The results from Severus are also
filtered with a set of germline structural variants VCF (from the Human Pangenome Reference Consortium) to remove false positives. The filtering was performed using `svpack`. `svpack` also provides simple annotation based on Ensembl GFF (v101), and the filtered VCF can be found in the folder `Severus_filtered_vcf`.

Severus also provides visualization of complex SV clusters. The HTML outputs of all clusters are compressed into a single zip at `Severus_cluster_plots`.

`AnnotSV` is used to further annotate the structural variants into a TSV file, (`AnnotatedSeverusSV`). The TSV file format is described in [AnnotSV README](https://github.com/lgmgeo/AnnotSV/blob/master/README.AnnotSV_latest.pdf). This provides more detailed annotation than `svpack`.

Lastly, to help with prioritizing variants relevant to cancer, the workflow also annotates the SVs with IntOGen Compendium of Cancer Genes (CCG) and produces a final set of structural variants in the `Annotated*SV_intogen` folder.

In the `SV_circos` folder, we generate a circos plot of all BND events that are at least 100 kbp apart. The circos plot is useful for visualizing large-scale rearrangements in the genome. Red color denotes fusion pairs that have been recorded in the Mitelman fusion database, whereas black color denotes fusions involving known genes in the GFF supplied to `svpack`. Grey color denotes fusions that are not in the Mitelman database or the GFF file.

## SNV/INDEL (Somatic and germline)

Somatic SNV/INDEL VCF can be found in the folder `small_variant_vcf`. We use DeepSomatic (trained on Revio dataset) or ClairS (default due to it being faster and trained on both Sequel II and Revio) to call both somatic SNV/INDELs. (See [input JSON parameters](step-by-step.md#input-json-parameters) 
for parameter to switch between them or to switch to Sequel IIe model for ClairS). 

DeepSomatic is currently computationally expensive (14-18 hours for 60X/30X tumor/normal) and requires using a separate caller for germline SNV/INDELs. The workflow implements Clair3 to call germline variants in both tumor and normal in addition to somatic variants. In addition, as DeepSomatic does not currently output the depth of coverage of the variants in normal, the VCF cannot be used as an input for Purple for purity and ploidy estimation. If DeepSomatic is used, Purple will run without somatic VCF which may affect estimation in some cases.

For both ClairS and DeepSomatic, The workflow can split the human genome into chunks (default 75 Mbp per chunk) and calls SNV/INDELs in parallel, then gathers the output into a single VCF. This allows the caller to scale to large genomes and large datasets by making use of multiple HPC nodes. Germline variants are called with Clair3 regardless of the somatic variant caller used.

For annotation, we use Ensembl VEP to annotate the VCF file (`small_variant_tsv_annotated` and `small_variant_vcf_annotated` folder). As with SV, the workflow also annotates the SNV/INDELs with IntOGen Compendium of Cancer Genes (CCG) and produce a final set of SNV/INDELs in the `small_variant_tsv_CCG` folder.

## Phasing

The workflow uses hiphase to phase both germline and somatic variants. As hiphase is a diploid phasing tool, somatic variants will be forced into parental haplotypes instead of the poly-clonal nature of the tumor. This is still useful for understanding the haplotype structure of the tumor. In the final VCF, `PS` tag represents the phase blocks. I.e. if two variants have the same `PS` tag, they belong to the same haplotype (within the phase block). In addition, we manually modified all the genotypes for called somatic variants to `0/1` (DeepSomatic by default assigns `1/1` to all somatic variants). This is done to allow HiPhase to phase the somatic variants.

## Homologous-recombination deficiency prediction

SNV, INDELs and SVs are supplied to CHORD for HRD prediction. The results can be found in the `chord_hrd_prediction` folder. We've tested internally with two datasets that have HRD and found that the results were accurate. However, the results have not been validated extensively with more samples and should be used with caution, especially at low tumor purity and coverage. E.g. at effective tumor coverage of 15X (30X tumor coverage with 50% tumor purity or 60X with 25% tumor purity), CHORD predicted HCC1395 to be HR-deficient.

## Differentially methylated region

CpG calls, at each loci in the human genome, are summarized using [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools). The bed file for the CpG calls are then used to call DMRs using DSS. The pileup bed file is found in the `pileup_(normal\/tumor)_bed` folder.

`DMR_annotated` contains 5 folders for each patient. Each folder represents differentially methylated region 
annotated differently. Annotation is done using [annotatr](https://bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html).

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

For each TSV file, the first 11 columns (`areaStat` being the last) are identical and 
are produced by [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html). 
- `meanMethyl1` is the mean methylation level in tumor.
- `meanMethyl2` is the mean methylation level in normal.
- `length` is the length of the DMR.
- `nCG` is the number of CpG sites in the DMR. By default the workflow requires at least 50 CpG sites in any DMR region.
- `areaStat` is the area statistic of the DMR. The larger the area statistic, the more significant the DMR is.

`annot.X` columns are produced by `annotatr` and all upper-case columns are extracted from the IntOGen Compendium of Cancer Genes TSV file.

## Somatic SNV/INDEL annotation

`small_variant_vcf_annotated` contains somatic variants annotated using Ensembl VEP. To obtain readable annotations, try using the `split-vep` tool from [`bcftools`](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html). The workflow
by default produces a TSV using `split-vep` in the `small_variant_tsv` folder.

## Prioritization

The workflow implements simple prioritization based on the IntOGen Compendium of Cancer Genes (CCG). The prioritization is done by
subsetting the DMR regions involving any gene in CCG with >50 CpG sites. For SV/SNV/INDELs, similar prioritization is done
by subsetting the variants involving any gene in CCG.

## Mutational signatures

Mutational signatures are determined using the `MutationalPattern` R [package](https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html).
There are two folders in the output:
- `mutsig_SNV` contains TSV files of the mutational signatures fitted to known COSMIC signatures. The TSVs are:
  - `mut_sigs.tsv` contains the contribution of the fitted signatures. The signature is fitted using the `fit_to_signatures_strict` function.
  - `reconstructed_sigs.tsv` contains the reconstructed trinucleotide context based on the fitted signatures.
  - `type_occurences.tsv` contains the frequency of each mutation type in the sample.
  - `mut_sigs_bootstrapped.tsv` contains the bootstrap results of the fitted signatures using the `fit_to_signatures_bootstrapped` function.
    The bootstrap results can be used to determine how stable the fitted signatures are.

## Purity and ploidy estimation

The workflow currently implements the HMFtools suite to estimate purity and ploidy based on Amber, Cobalt and Purple for tumor/normal workflow. 
For Cobalt, due to noises in long-reads based read-depth and B-allele frequency, we set PCF gamma to 1000 to allow for 
better segmentation. See [GitHub issue](https://github.com/hartwigmedical/hmftools/issues/485) for discussion. Practically, this means
that the CNV calls may miss smaller focal events, but those should be picked up by the SV caller.

The purity and ploidy estimates were found to be reasonably robust in our experience with HCC1395 and COLO829, 
but should be used with caution. Purity and ploidy estimates can be found in the `*.purity.tsv` file in `Purple_outputs` folder.

Purity and ploidy are also estimated by Wakhan, see section below on copy number variation.

## Copy number variation

CNVKit is used to segment copy numbers from the matched tumor/normal BAM files. To optimize for long-reads, we set bin size to 10 kbp and found
it to be optimal based on COLO829. The workflow also uses purity and ploidy estimates from HMFtools in combination with ClairS heterozygous SNVs
to estimate allele-specific major and minor copy numbers in the `cnvkit_cns_with_major_minor_CN` folder. 

A downside of CNVKit is that the recall mode requires integer ploidy, which can fail in cases where there are subclonal CNV.
Purple also calls allele-specific copy numbers and is able to account for non-integer ploidy (subclonal CNV). The CNV segments from Purple can be found in the `Purple_outputs` folder and has the suffix `tumor.purple.cnv.somatic.tsv`. We visualize the results from Purple in the final report.

Lastly, we also include [Wakhan](https://github.com/KolmogorovLab/Wakhan) as an alternative CNV caller for both tumor/normal and tumor-only workflow. Wakhan calls haplotype-specific CNV and also estimates the purity and ploidy of the sample. The results can be found in the `wakhan_cnv` folder. We only choose the best purity and ploidy solutions in the final results folder, but the full results including alternative solutions picked by Wakhan can be found by extracting the `tar.gz` file. As of version 0.9, Wakhan is still in active development so we report `Purple` results by default for tumor/normal workflow. For tumor-only workflow, `Wakhan` is used in the final report as Purple has not been tested for tumor-only workflow.

## Report html

The workflow produces a report html file that summarizes the results from the workflow. The report html file can be found in the `report` folder. The report also contains variants filtered with a set of criteria documented in the HTML file to help with the interpretation of the samples. However, the filtering criteria are not meant to be used as a hard rule and should be used with caution. They may not be optimal for all datasets and we encourate users to analyse the pipeline results in more detail.
