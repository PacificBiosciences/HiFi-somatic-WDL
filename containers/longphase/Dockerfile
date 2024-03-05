FROM continuumio/miniconda3:23.10.0-1

RUN apt-get update && \
    apt-get install -y git g++ gcc autoconf make zlib1g-dev libbz2-dev liblzma-dev procps git libarchive-dev && \
    rm -rf /var/lib/apt/lists/* 

RUN conda install -n base -c conda-forge mamba libarchive
RUN mamba create -n tools
RUN mamba install -c bioconda -c conda-forge -n tools csvtk==0.27.2 samtools==1.17 bcftools==1.17 bedtools==2.31.0 datamash==1.8 seqkit==2.5.1
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate tools" >> ~/.bashrc
RUN conda clean --all -y

WORKDIR /opt/longphase
RUN git clone https://github.com/twolinin/longphase.git /opt/longphase && \
    git checkout 6c5a9c798a9daa0f77d14bff35f20c4424bfb6ec && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase/.git

ENV PATH=/opt/conda/envs/tools/bin:${PATH}:${HOME}/bin:/opt/longphase

CMD ["longphase", "phase", "--help"]
