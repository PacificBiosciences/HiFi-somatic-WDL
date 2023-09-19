FROM continuumio/miniconda3

WORKDIR /app

RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update
RUN apt-get install -y procps git
RUN conda install -n base -c conda-forge mamba libarchive
RUN mamba create -n severus
RUN mamba install -c bioconda -c conda-forge -n severus numpy==1.25.1 networkx==3.1 pydot==1.4.2 graphviz==8.1.0 pysam==0.21.0 samtools==1.17 bcftools==1.17 htslib==1.17
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate severus" >> ~/.bashrc
ENV PATH /opt/conda/envs/severus/bin:$PATH
RUN git clone https://github.com/KolmogorovLab/Severus.git
RUN pip cache remove *
RUN conda clean --all
