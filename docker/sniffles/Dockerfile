FROM ubuntu:22.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    git \
    wget \
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

WORKDIR /htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar jxvf htslib-1.19.tar.bz2  
WORKDIR /htslib/htslib-1.19
RUN ./configure && \
    make && \
    make install 
    
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
