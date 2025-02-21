# Builder stage: Install dependencies and create the conda environment
FROM continuumio/miniconda3:23.10.0-1 AS builder

WORKDIR /app

# Install system dependencies and minimal tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps git libarchive-dev && \
    rm -rf /var/lib/apt/lists/*

# Install mamba and build the 'tools' environment
RUN conda install -n base -c conda-forge mamba libarchive && \
    mamba create -n tools -y && \
    mamba install -n tools -c bioconda -c conda-forge \
      csvtk==0.27.2 \
      samtools==1.17 \
      bcftools==1.17 \
      bedtools==2.31.0 \
      datamash==1.8 \
      seqkit==2.5.1 \
      pycirclize \
      pandas && \
    conda clean --all --yes

# Copy application files into the builder stage
COPY Compendium_Cancer_Genes_2024-9-30.tsv /app/Compendium_Cancer_Genes.tsv
COPY hg38.bed hg38_cytoband.tsv plot_circos.py MCGENE.TXT.DATA /app/

# Final stage: Use a fresh image and copy only necessary artifacts
FROM continuumio/miniconda3:23.10.0-1

WORKDIR /app

# Copy the built conda environment and application files from the builder
COPY --from=builder /opt/conda/envs/tools /opt/conda/envs/tools
COPY --from=builder /app/ /app/

# Set up environment to use the 'tools' environment
ENV PATH /opt/conda/envs/tools/bin:$PATH

CMD ["bash"]