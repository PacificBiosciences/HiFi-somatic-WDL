FROM bioconductor/bioconductor_docker

WORKDIR /app

RUN apt-get update 
RUN apt-get install -y default-jre python3 python3-pip wget
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*
RUN R -e 'BiocManager::install(c("tidyverse", "VariantAnnotation", "copynumber"))'

RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/amber-v3.9.1/amber-3.9.1.jar'
RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.15.1/cobalt_v1.15.1.jar'
RUN wget 'https://github.com/hartwigmedical/hmftools/releases/download/purple-v3.9/purple_v3.9.jar'
RUN wget 'https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz'

RUN mkdir -p /app/resources
RUN tar xfz hmf_dna_pipeline_resources.38_v5.33.tar.gz -C /app/resources
RUN rm -f hmf_dna_pipeline_resources.38_v5.33.tar.gz