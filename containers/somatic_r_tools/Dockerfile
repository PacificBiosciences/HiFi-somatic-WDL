FROM bioconductor/bioconductor_docker:RELEASE_3_17

WORKDIR /app

RUN apt-get update 
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN R -e 'BiocManager::install(c("DSS", "data.table", "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "cowplot", "ggplot2", "GenomicRanges", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "annotatr"))'
COPY DSS_tumor_normal.R .
COPY mutational_pattern.R .
COPY annotatr_dmr.R .
