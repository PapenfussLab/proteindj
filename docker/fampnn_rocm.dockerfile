FROM quay.io/sarahbeecroft9/rocm-mpich-base:rocm7.2.3-mpich3.4.3-ubuntu24.04
# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update --quiet \
    && apt-get install --no-install-recommends --yes --quiet \
        build-essential \
        git \
        wget \
        python-dev-is-python3 \
        python3-venv \
        python3-pip \
        build-essential \
        software-properties-common \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove --yes \
    && apt-get clean

# Install miniforge (mamba)
RUN set -eux ; \
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" ; \
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3 -s ; \
    rm -rf ./Miniforge3-*
ENV PATH /opt/miniforge3/bin:$PATH

RUN mamba install -c conda-forge -c bioconda python=3.10 \
    && mamba clean -afy
# clone repository
RUN git clone https://github.com/PapenfussLab/fampnn.git /app/fampnn
RUN pip install torch==2.12.1 torchvision==0.27.1 --index-url https://download.pytorch.org/whl/rocm7.2

RUN pip install scipy \
    certifi \
    torch-geometric \
    gemmi \
    biopython \
    tqdm \
    natsort \
    omegaconf \
    rdkit \
    einops \
    hydra-core \
    torchtyping \
    dm-tree \
    timm \
    jaxtyping \
    joblib \
    pandas

RUN cd /app/fampnn \
    && pip install -e .

ENV PYTHONPATH="/app/fampnn:$PYTHONPATH"

RUN git clone https://github.com/SarahBeecroft/proteindj.git /tmp/proteindj \
    && mv /tmp/proteindj/scripts /scripts \
    && rm -rf /tmp

COPY fampnn_rocm.dockerfile /opt/docker-recipes/
LABEL org.opencontainers.image.authors="Sarah Beecroft <sarah.beecroft@csiro.au>" \
      org.opencontainers.image.description="Fampnn ROCm7.2.3 support on Ubuntu 24.04 for AMD GPUs" \
      Python_version="3.10"
