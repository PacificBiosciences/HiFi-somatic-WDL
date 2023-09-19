FROM bioconductor/bioconductor_docker

WORKDIR /app

RUN apt-get update 
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

RUN R -e 'BiocManager::install(c("DSS", "data.table"))'
COPY DSS_tumor_normal.R .
