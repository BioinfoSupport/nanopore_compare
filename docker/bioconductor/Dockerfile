FROM docker://bioconductor/bioconductor_docker:3.18

RUN Rscript -e 'BiocManager::install(c("data.table", "R.utils", "GenomicAlignments", "VariantAnnotation"))'

# RUN	hts_ver=1.19; \
# 	for p in samtools bcftools htslib; do \
# 		curl -O -L https://github.com/samtools/$p/releases/download/${hts_ver}/$p-${hts_ver}.tar.bz2; \
# 		tar xvjf $p-${hts_ver}.tar.bz2; \
# 	        cd $p-${hts_ver}; \
# 		./configure; \
# 		make -j; \
# 		make install; \
# 		cd ..; \
# 		rm -rf $p-${hts_ver}*; \
# 	done

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

################## METADATA ######################
LABEL maintainer="fedor.bezrukov@gmail.com"
