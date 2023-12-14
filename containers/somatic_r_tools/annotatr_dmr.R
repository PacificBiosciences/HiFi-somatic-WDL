library(annotatr)
library(GenomicRanges)
library(data.table)

# Read in argument. Just the DMRs output from DSS
args <- commandArgs(trailingOnly = TRUE)
# Fail if not enough arguments
if (length(args) < 3) {
    stop("Not enough arguments. Need DMR file, patient name and number of threads.")
}

dmr_file <- args[1]
patient_name <- args[2]
threads <- as.integer(args[3])

dmrs <- fread(dmr_file, nThread = threads)
# Convert to Granges
dmrs <- GRanges(dmrs)

ann_to_build <- c(
    "hg38_genes_promoters", "hg38_genes_1to5kb",
    "hg38_genes_5UTRs", "hg38_genes_exons",
    "hg38_genes_introns", "hg38_genes_3UTRs"
)

# Build annotations for hg38
built_annotations <- build_annotations(genome = "hg38", annotations = ann_to_build)

# Annotate DMRs
dmrs <- annotate_regions(dmrs, built_annotations, ignore.strand = TRUE)

# Summarize annotation type
dm_annsum <- summarize_annotations(
    annotated_regions = dmrs,
    quiet = TRUE
)

colnames(dm_annsum) <- c("Annotation_Type", "Count")

# Write a table of the annotation summary
fwrite(
    dm_annsum,
    paste0(patient_name, "_dmr_annotation_summary.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA",
    append = FALSE, 
    compress = "gzip",
    nThread = threads
)

# Split by annotation type and write each type to a file
for (ann_type in ann_to_build) {
    # Extract from dmrs granges annot$type column and save
    # to a TSV file
    # Need to unroll Granges at annot column
    ann_type_dmrs <- as.data.table(dmrs[dmrs$annot$type == ann_type])
    # Write to a file
    fwrite(
        ann_type_dmrs,
        paste0(patient_name, "_", ann_type, "_dmrs.tsv.gz"),
        sep = "\t",
        quote = FALSE,
        na = "NA",
        append = FALSE,
        compress = "gzip",
        nThread = threads
    )
}

sessionInfo()
version