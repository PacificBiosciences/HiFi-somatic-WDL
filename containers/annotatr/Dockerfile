FROM bioconductor/bioconductor_docker

WORKDIR /app

RUN apt-get update 
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN R -e 'BiocManager::install(c("GenomicRanges", "data.table", "annotatr", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"))'

COPY annotatr_dmr.R /app/annotatr_dmr.R