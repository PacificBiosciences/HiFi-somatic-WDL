FROM bioconductor/bioconductor_docker:RELEASE_3_14

WORKDIR /app

RUN apt-get update 
RUN apt-get install -y default-jre python3 python3-pip wget
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*
RUN R -e 'BiocManager::install(c("tidyverse", "VariantAnnotation", "copynumber"))'

RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/amber-v3.9.1/amber-3.9.1.jar'
RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.14.1/cobalt-1.14.1.jar'
RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/purple-v3.9/purple_v3.9.jar'
