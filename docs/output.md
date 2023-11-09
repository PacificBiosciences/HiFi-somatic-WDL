# Output from workflow

- [Output from workflow](#output-from-workflow)
  - [QC results](#qc-results)
  - [Structural variants](#structural-variants)
  - [SNV/INDEL (Somatic and germline)](#snvindel-somatic-and-germline)
  - [Differentially methylated region](#differentially-methylated-region)
  - [Somatic SNV/INDEL annotation](#somatic-snvindel-annotation)
  - [Prioritization](#prioritization)
  - [Mutational signatures](#mutational-signatures)
  - [Purity and ploidy estimation](#purity-and-ploidy-estimation)
  - [Copy number variation](#copy-number-variation)

## QC results

Depth of coverage for both tumor and normal can be found in the `mosdepth_normal_summary` and `mosdepth_tumor_summary` folders. 
The `mosdepth` folder contains the depth of coverage for each chromosome and overall coverage (`total` at the end of the file).

`overall_(tumor\|normal)_alignment_stats` contains the alignment statistics for both tumor and normal. The statistics are produced by
using `seqkit bam` command. Per-alignment stats are also generated in the `per_alignment_(tumor\|normal)_stats` folders.

## Structural variants

The pipeline incorporates two structural variants callers: [Severus](https://github.com/KolmogorovLab/Severus) and [Sniffles](https://github.com/fritzsedlazeck/Sniffles). 
For Sniffles, we found version 2.0.7 to have higher recall rate for the COLO829 benchmark dataset. The workflow for Sniffles is as follow:

- Generate SNF for tumor and normal.
- Joint-call tumor and normal SNF without any QC. 
- The VCF is then filtered using slivar to >0 coverage in normal, 0 reference reads and >= 2 variant reads. 
- Filtered the resultant VCF with a Sniffles joint-called cohort of 110+ HPRC samples. 

For Severus, the only additional parameter modified is `--min-support 2` to maximize the sensitivity of the caller. The results from Severus are also
filtered with the HPRC control VCF to remove false positives. For both callers, the filtering was done with Truvari using the following command:

``` bash
truvari bench -p 0 -s 0 -S 0 --sizemax 100000000 --dup-to-ins
```

The `fp.vcf.gz` VCF file from the `bench` step is taken as the final VCF file. This is why the final SV VCF files contain Truvari annotations in the INFO field (folder `Severus_filterHPRC_vcf` and `sniffles_somatic_vcf_filterHPRC`).

We then use AnnotSV to annotate the structural variants into a TSV file (`AnnotatedSeverusSV` and `AnnotatedSnifflesSV` folders). The TSV file format is described in [AnnotSV README](https://github.com/lgmgeo/AnnotSV/blob/master/README.AnnotSV_latest.pdf).
To help prioritizing variants relevant to cancer, the workflow also annotates the SVs with IntOGen Compendium of Cancer Genes (CCG) and produce a final set of SV in the `Annotated*SV_intogen` folder

## SNV/INDEL (Somatic and germline)

We use ClairS to call both somatic and germline SNV/INDELs (Default Revio model, please see [input JSON parameters](step-by-step.md#input-json-parameters) 
for parameter to switch to Sequel II). The workflow splits the human genome into chunks (default 50 Mbp) and calls SNV/INDELs in parallel, then collect
the output into a single VCF (`small_variant_vcf` folder). Finally, we use Ensembl VEP to annotate the VCF file (`small_variant_tsv_annotated` folder). Similar to SV, the workflow also annotates the SNV/INDELs with IntOGen Compendium of Cancer Genes (CCG) and produce a final set of SNV/INDELs in the `small_variant_tsv_CCG` folder.

## Differentially methylated region

Firstly, CpG call at each loci in the human genome is summarized using [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools). The bed file for the CpG call is then used to call DMRs using DSS. You
may find the pileup bed file in the `pileup_(normal\/tumor)_bed` folder.

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
- `meanMethyl1` refers to the mean methylation level in tumor.
- `meanMethyl2` refers to the mean methylation level in normal.
- `length` refers to the length of the DMR.
- `nCG` refers to then number of CpG sites in the DMR. By default the workflow requires at least 50 CpG sites in any DMR region.
- `areaStat` refers to the area statistic of the DMR. The larger the area statistic, the more significant the DMR is.

`annot.X` columns are produced by `annotatr` and all upper-case columns are extracted from IntOGen Compendium of Cancer Genes TSV file.

## Somatic SNV/INDEL annotation

`small_variant_vcf_annotated` contains somatic variants annotated using Ensembl VEP. In order to obtain readable annotation 
you may try using the `split-vep` tool from [`bcftools`](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html). The workflow
by default produces a TSV using `split-vep` in the `small_variant_tsv` folder.

## Prioritization

The workflow implements simple prioritization based on IntOGen Compendium of Cancer Genes (CCG). The prioritization is done by
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

The workflow currently implements the HMFtools suite to estimate purity and ploidy based on Amber, Cobalt and Purple. However, 
Cobalt does not currently count the reads in long-reads correctly so the read-depth segmentation is not accurate. Nevertheless, the purity
and ploidy estimates were found to be reasonably robust in our experience, but should be used with caution. Purity and ploidy estimates can be found
in the `*.purity.tsv` file in `Purple_outputs` folder.

## Copy number variation

CNVKit is used to segment copy number from the matched tumor/normal BAM files. To optimize for long-reads, we set bin size to 10 kbp and found
it to be optimal based on COLO829. The workflow also uses purity and ploidy estimates from HMFtools in combination with ClairS heterozygous SNVs
to estimate major and minor copy numbers in `cnvkit_cns_with_major_minor_CN` folder.