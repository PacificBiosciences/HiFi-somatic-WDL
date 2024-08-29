FROM bioconductor/bioconductor_docker:RELEASE_3_14

WORKDIR /app

RUN apt-get update && apt-get install -y default-jre python3 python3-pip wget maven git circos && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN R -e 'BiocManager::install(c("tidyverse", "VariantAnnotation", "copynumber"))'
RUN R -e 'BiocManager::install(c("cowplot"))'

RUN wget https://github.com/hartwigmedical/hmftools/archive/refs/tags/amber-v4.0.tar.gz
RUN wget https://github.com/hartwigmedical/hmftools/archive/refs/tags/cobalt-v1.16.tar.gz
RUN wget https://github.com/hartwigmedical/hmftools/archive/refs/tags/purple-v4.0.tar.gz

# Copy jar file. These were built from hash 25863b0367f58dba6dbde64223ffbd9fbad96f28
# with maven 3.6.3 and jdk 11
COPY amber-4.0-jar-with-dependencies.jar amber.jar
COPY amber-4.0-jar-with-dependencies.gamma1000.jar amber.gamma1000.jar
COPY cobalt-1.16.0-jar-with-dependencies.jar cobalt.jar
COPY purple-4.0-jar-with-dependencies.jar purple.jar
COPY process_clairS_VCF.sh process_clairS_VCF.sh

# Circos
COPY circos-0.69-9.tgz circos-0.69-9.tgz
RUN tar -xzf circos-0.69-9.tgz && mkdir /usr/bin/etc/ && cp -r circos-0.69-9/data /usr/bin/etc/ && ln -s /etc/circos/* /usr/bin/etc/ && cp -r circos-0.69-9/error /usr/bin/
