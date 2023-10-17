FROM continuumio/miniconda3

WORKDIR /app

RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update
RUN apt-get install -y procps git libarchive-dev
RUN conda install -n base -c conda-forge mamba libarchive
RUN mamba create -n tools
RUN mamba install -c bioconda -c conda-forge -n tools csvtk==0.27.2 samtools==1.17 bcftools==1.17 bedtools==2.31.0 datamash==1.8 seqkit==2.5.1
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate tools" >> ~/.bashrc
ENV PATH /opt/conda/envs/tools/bin:$PATH
RUN pip cache remove *
RUN conda clean --all
COPY Compendium_Cancer_Genes.tsv /app/Compendium_Cancer_Genes.tsv
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*