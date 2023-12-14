library(DSS)
library(bsseq)
library(data.table)

# Read in arguments. First for tumor bed file, then normal bed
# file, then output file name
args <- commandArgs(trailingOnly = TRUE)
# Fail if not enough arguments
if (length(args) < 3) {
    stop("Not enough arguments. Need tumor bed file, normal bed file, and output file name.")
}
tumor_file <- args[1]
normal_file <- args[2]
output_file <- args[3]
threads <- args[4]
# Hardcode this for now
threads <- 8

# Read in count table
tumor <- fread(tumor_file, header = FALSE, sep = "\t")
normal <- fread(normal_file, header = FALSE, sep = "\t")
colnames(tumor) <- c("chr", "pos", "N", "X")
colnames(normal) <- c("chr", "pos", "N", "X")

bs_obj <- makeBSseqData(
    list(
        tumor,
        normal
    ),
    sampleNames = c("tumor", "normal")
)

rm(list = c("tumor", "normal"))

dmlTest.sm <- DMLtest(
    bs_obj,
    group1 = c("tumor"), group2 = c("normal"),
    smoothing = TRUE, smoothing.span = 500,
    ncores = threads
)

rm(list = c("bs_obj"))

dmrs <- callDMR(
    dmlTest.sm,
    p.threshold = 0.001,
    delta = 0.25,
    minCG = 5,
    minlen = 100
)

fwrite(
    dmrs,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    append = FALSE, 
    nThread = threads
)

## Loop through chromosomes
# # Split into chromosome
# chromosomes_tumor <- unique(tumor$chr)
# chromosomes_normal <- unique(normal$chr)
# # Get common chromosome
# chromosomes <- intersect(chromosomes_tumor, chromosomes_normal)

# count <- 0
# for (chrom in chromosomes) {
#     print(paste0("Processing chr: ", chrom))

#     bs_obj <- makeBSseqData(
#         list(
#             tumor[chr == chrom],
#             normal[chr == chrom]
#         ),
#         sampleNames = c("tumor", "normal")
#     )

#     dmlTest.sm <- DMLtest(
#         bs_obj,
#         group1 = c("tumor"), group2 = c("normal"),
#         smoothing = TRUE, smoothing.span = 500,
#         ncores = 8
#     )

#     # dmls = callDML(
#     #     dmlTest, p.threshold=0.001,
#     #     delta = 0.25
#     #     )

#     dmrs <- callDMR(
#         dmlTest.sm,
#         p.threshold = 0.001,
#         delta = 0.25,
#         minCG = 5,
#         minlen = 100
#     )

#     # Save dmrs with data.table
#     if (count == 0 && !is.null(dim(dmrs))) {
#         fwrite(
#             dmrs,
#             output_file,
#             sep = "\t",
#             quote = FALSE,
#             na = "NA",
#             append = FALSE
#         )
#     } else if (count > 0 && !is.null(dim(dmrs))) {
#         fwrite(
#             dmrs,
#             output_file,
#             sep = "\t",
#             quote = FALSE,
#             na = "NA",
#             append = TRUE
#         )
#     } else {
#         print(paste0("No DMRs found for chr: ", chrom))
#     }

#     count <- count + 1
# }

sessionInfo()
version