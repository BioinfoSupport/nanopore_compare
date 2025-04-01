FROM ubuntu:22.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    git \
    curl \
    bzip2 \
    make \
    cmake \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3.10 \
    python3.10-venv \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

LABEL samtools="1.21"
ARG hts_ver=1.21
WORKDIR /htslib
RUN curl -L https://github.com/samtools/htslib/releases/download/$hts_ver/htslib-$hts_ver.tar.bz2 | tar -xjf -
WORKDIR /htslib/htslib-${hts_ver}
RUN ./configure && \
    make -j && \
    make install 
WORKDIR /htslib
RUN curl -L https://github.com/samtools/bcftools/releases/download/$hts_ver/bcftools-$hts_ver.tar.bz2 | tar -xjf -
WORKDIR /htslib/bcftools-${hts_ver}
RUN ./configure && \
    make -j && \
    make install 
WORKDIR /
RUN rm -rf htslib

# RUN pip3 install sniffles
RUN python3 -m venv /venv \
  && pip install --no-cache-dir pysam==0.22.0 \
  && pip install --no-cache-dir sniffles==2.2

################## METADATA ######################
LABEL maintainer="fedor.bezrukov@gmail.com"
LABEL version=2.2
LABEL pysam=0.22.0
LABEL htslib=1.19

CMD sniffles
