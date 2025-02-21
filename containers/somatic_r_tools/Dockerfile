FROM bioconductor/bioconductor_docker:RELEASE_3_18 AS builder

WORKDIR /app

RUN R -e 'BiocManager::install(c("DSS", "data.table", "MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38", "cowplot", "ggplot2", "GenomicRanges", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "annotatr", "DT", "quarto", "qualpalr", "dplyr"))' &&\
    R -e 'install.packages("tinytable", repos = "https://vincentarelbundock.r-universe.dev")'

# Stage 2: Final image
FROM bioconductor/bioconductor_docker:RELEASE_3_18

WORKDIR /app

COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.56/quarto-1.5.56-linux-amd64.deb &&\
    dpkg -i quarto-1.5.56-linux-amd64.deb &&\
    rm -rf quarto-1.5.56-linux-amd64.deb /usr/lib/rstudio-server &&\
    find /usr/local/lib/R/site-library -name '*.so' -exec strip {} \; &&\
    mkdir -p /app/references

COPY DSS_tumor_normal.R mutational_pattern.R annotatr_dmr.R visualize_hifisomatic.qmd visualize_hifisomatic_tumor_only.qmd ./
COPY cytoBand_hg38.txt.gz copy_number_functions.R common_func.R /app/references/

