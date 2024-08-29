FROM bioconductor/bioconductor_docker:RELEASE_3_18

WORKDIR /app

RUN apt-get update &&\
    apt-get clean &&\
    rm -rf /var/lib/apt/lists/* &&\
    mkdir -p /app/references

RUN R -e 'BiocManager::install(c("DSS", "data.table", "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "cowplot", "ggplot2", "GenomicRanges", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "annotatr", "DT", "quarto", "qualpalr", "dplyr"))' &&\
    R -e 'install.packages("tinytable", repos = "https://vincentarelbundock.r-universe.dev")'

# R -e 'remotes::install_github("vincentarelbundock/tinytable@d0a2d1f38dee5907750781c0a50967a7507e8b22")'

COPY DSS_tumor_normal.R mutational_pattern.R annotatr_dmr.R visualize_hifisomatic.qmd ./
COPY cytoBand_hg38.txt.gz copy_number_functions.R /app/references/

RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.56/quarto-1.5.56-linux-amd64.deb &&\
    dpkg -i quarto-1.5.56-linux-amd64.deb &&\
    rm quarto-1.5.56-linux-amd64.deb
