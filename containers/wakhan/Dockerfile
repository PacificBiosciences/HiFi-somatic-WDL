# Optimized single-stage build for Wakhan
FROM mambaorg/micromamba:1.5.3

USER root

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_ENV_NAME=wakhan

# Install git and minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
        procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Switch back to mambauser
USER mambauser

# Create /opt/wakhan directory and change ownership
USER root
RUN mkdir -p /opt/wakhan && chown -R mambauser:mambauser /opt/wakhan
USER mambauser

# Clone Wakhan repository at specific commit
RUN git clone --depth 1 https://github.com/KolmogorovLab/Wakhan.git /tmp/wakhan && \
    cd /tmp/wakhan && \
    git fetch --depth 1 origin 2c909df59337f6c52c8437db5e934d2aee06ceec && \
    git checkout 2c909df59337f6c52c8437db5e934d2aee06ceec && \
    cp -r /tmp/wakhan/* /opt/wakhan/

# Create optimized environment with specific versions and install Wakhan
RUN echo "name: wakhan" > /tmp/environment.yml && \
    echo "channels:" >> /tmp/environment.yml && \
    echo "  - conda-forge" >> /tmp/environment.yml && \
    echo "  - bioconda" >> /tmp/environment.yml && \
    echo "dependencies:" >> /tmp/environment.yml && \
    echo "  - python=3.8" >> /tmp/environment.yml && \
    echo "  - samtools>=1.14" >> /tmp/environment.yml && \
    echo "  - bcftools>=1.14" >> /tmp/environment.yml && \
    echo "  - pip" >> /tmp/environment.yml && \
    echo "  - pip:" >> /tmp/environment.yml && \
    echo "    - pysam" >> /tmp/environment.yml && \
    echo "    - pyfaidx" >> /tmp/environment.yml && \
    echo "    - numpy" >> /tmp/environment.yml && \
    echo "    - pandas" >> /tmp/environment.yml && \
    echo "    - plotly" >> /tmp/environment.yml && \
    echo "    - scikit-learn==1.2.0" >> /tmp/environment.yml && \
    echo "    - scipy==1.9.2" >> /tmp/environment.yml && \
    echo "    - ruptures" >> /tmp/environment.yml && \
    echo "    - vcf_parser" >> /tmp/environment.yml && \
    echo "    - kaleido" >> /tmp/environment.yml && \
    micromamba create -f /tmp/environment.yml -n $CONDA_ENV_NAME -y && \
    micromamba clean --all --yes && \
    eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate $CONDA_ENV_NAME && \
    cd /opt/wakhan && \
    pip install -e . --no-cache-dir && \
    rm -rf /tmp/wakhan /tmp/environment.yml

# Set working directory
WORKDIR /opt/wakhan

# Set environment variables for runtime
ENV PATH="/opt/conda/envs/$CONDA_ENV_NAME/bin:$PATH"

# Activate conda environment by default
RUN echo "micromamba activate $CONDA_ENV_NAME" >> ~/.bashrc

# Set entry point to use wakhan command
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/wakhan/wakhan.py"]

# Default command shows help
CMD ["--help"]

# Add labels for metadata
LABEL maintainer="Wakhan Docker Image"
LABEL version="commit 2c909df59337f6c52c8437db5e934d2aee06ceec"
LABEL description="Wakhan - A tool for analyzing haplotype-specific chromosome-scale somatic copy number aberrations"
LABEL repository="https://github.com/KolmogorovLab/Wakhan"
