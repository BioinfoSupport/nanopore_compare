# Building paf2chain stage
FROM docker://rust:1.76.0 as builder

WORKDIR /usr/src/paf2chain
COPY paf2chain .
RUN cargo install --path .


# Main stage
FROM docker://bioconductor/bioconductor_docker:3.18

RUN Rscript -e 'BiocManager::install(c("data.table", "R.utils", "GenomicAlignments", "VariantAnnotation"))'

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    bcftools \
    bzip2 \
    cmake \
    gcc \
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
    python3.10-venv \
    tabix \
    wget \
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
## Do this from workflow cache ?
# && \
#     cd /usr/local/lib64/python3.10/site-packages/medaka/data/ && \
#     curl -L -O https://github.com/nanoporetech/medaka/raw/master/medaka/data/r1041_e82_400bps_sup_g615_model.tar.gz

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

# Copy paf2chain from the build stage
COPY --from=builder /usr/local/cargo/bin/paf2chain /usr/local/bin/paf2chain

RUN pip3 install CrossMap && \
    rm -rf /root/.cache

WORKDIR /

	## Cleanup
	# Not needed in the image
	#dnf erase -y cmake gcc-g++ make zlib-devel git python3-devel ncurses-devel bzip2-devel xz-devel
#	dnf clean all
