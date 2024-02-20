### Building paf2chain stage
FROM docker://rust:1.76.0 as rust_builder

WORKDIR /usr/src/paf2chain
COPY paf2chain .
RUN cargo install --path .

### Build the rest: Flye and samtools/htslib
FROM docker://ubuntu:22.04 as builder

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    bzip2 \
    cmake \
    curl \
    gcc \
    g++ \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    make \
    ncurses-dev \
    python-is-python3 \
    python3-pip \
    python3.10 \
    python3.10-dev \
    python3.10-venv \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /
RUN git clone -b 2.9.3 --depth 1 https://github.com/fenderglass/Flye
WORKDIR	/Flye
RUN make && \
    python3 setup.py install && \
    ln -s flye-minimap2 /usr/local/bin/minimap2 && \
    curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
    cp k8-0.2.4/k8-`uname -s` /usr/local/bin/k8 && \
    cp lib/minimap2/misc/*.js /usr/local/bin
WORKDIR /
RUN rm -rf Flye


RUN pip3 install "medaka==1.11.2" && \
    rm -rf /root/.cache

ARG hts_ver=1.19
WORKDIR /htslib
RUN curl -L https://github.com/samtools/htslib/releases/download/$hts_ver/htslib-$hts_ver.tar.bz2 | tar -xjf -
WORKDIR /htslib/htslib-${hts_ver}
RUN ./configure && \
    make -j && \
    make install 
WORKDIR /htslib
RUN curl -L https://github.com/samtools/samtools/releases/download/$hts_ver/samtools-$hts_ver.tar.bz2 | tar -xjf -
WORKDIR /htslib/samtools-${hts_ver}
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



### Main stage
FROM docker://ubuntu:22.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    bzip2 \
    curl \
    python-is-python3 \
    python3-pip \
    python3.10 \
    python3.10-venv \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Copy paf2chain from the build stage
COPY --from=rust_builder /usr/local/cargo/bin/* /usr/local/bin/
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/lib/python3.10/dist-packages /usr/local/lib/python3.10/dist-packages

RUN pip3 install CrossMap && \
    rm -rf /root/.cache

WORKDIR /
