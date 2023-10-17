FROM bioconductor/bioconductor_docker

WORKDIR /app

RUN apt-get update 
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN R -e 'BiocManager::install(c("DSS", "data.table", "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "cowplot", "ggplot2"))'
COPY DSS_tumor_normal.R .
COPY mutational_pattern.R .
