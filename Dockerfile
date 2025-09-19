FROM continuumio/miniconda3

LABEL maintainer="Samuel Bunga"
WORKDIR /app

# Install build tools and dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        zlib1g-dev \
        cmake \
        vim \
        git \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages into base
RUN conda install -y -c bioconda -c conda-forge \
    python=3.12 \
    biopython=1.85 \
    bedtools=2.31.1 \
    blast=2.16.0 \
    samtools=1.21 \
    matplotlib=3.10.6 \
    pandas=2.3.2 \
    seaborn=0.13.2 \
    multiqc=1.31 \
    bbmap=39.26 \
    openjdk=17 \
    && conda clean -afy

ENV JAVA_HOME=/opt/conda

# Install Nextflow
RUN wget -qO /usr/local/bin/nextflow https://github.com/nextflow-io/nextflow/releases/download/v25.04.2/nextflow \
    && chmod +x /usr/local/bin/nextflow

# Install SortMeRNA 4.3.6 via binary installer
RUN mkdir -p /opt/sortmerna && \
    chmod 777 /opt/sortmerna && \
    wget https://github.com/biocore/sortmerna/releases/download/v4.3.6/sortmerna-4.3.6-Linux.sh && \
    bash sortmerna-4.3.6-Linux.sh --skip-license --exclude-subdir --prefix=/opt/sortmerna && \
    rm sortmerna-4.3.6-Linux.sh

# Add SortMeRNA to PATH
ENV PATH="/opt/sortmerna/bin:$PATH"

# Confirm installed tools
RUN echo "Installed tools:" && \
    python --version && \
    bedtools --version && \
    blastn -version && \
    samtools --version && \
    sortmerna --version

# Copy the download script into the container
COPY scripts/download_sortmerna_db.sh /opt/download_sortmerna_db.sh

# Make it executable and run it
RUN chmod +x /opt/download_sortmerna_db.sh && \
    /opt/download_sortmerna_db.sh && \
    rm /opt/download_sortmerna_db.sh

# ENV variables
RUN mkdir -p '/app/logs'
ENV NXF_VER='25.04.2'
ENV NXF_LOG_FILE='/app/logs/nextflow.logs'
ENV NXF_CACHE_DIR='/app/nextflow_files/'

# Set working directory
WORKDIR /app
COPY ./bin /app/bin
COPY ./subworkflows /app/subworkflows
COPY ./main.nf /app/main.nf
COPY ./nextflow.config /app/nextflow.config
COPY ./config /app/config

SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]
