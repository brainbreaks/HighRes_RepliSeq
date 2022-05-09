FROM rocker/tidyverse:4.1.0
LABEL version="1"
LABEL software="Repli-seq"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

ARG http_proxy
ARG https_proxy
ENV DESTINATION=/bin
ENV MAKEFLAGS='-j 32'

# Install required standard libraries
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get update
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get install \
    locales locales-all python2 \
    make gcc gfortran g++ perl libcurl4-openssl-dev \
    ca-certificates ghostscript \
    zip unzip wget mime-support \
    picard-tools trim-galore \
    libopenblas-base libopenblas-dev liblapack-dev liblapack-dev \
    libncurses6 libncursesw6 libncurses-dev \
    zlib1g-dev bzip2 libbz2-dev liblzma-dev bowtie2 -y && \
    apt-get clean && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales
RUN ln -s /usr/bin/python2 /usr/bin/python
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" install2.r --error --deps TRUE optparse dplyr reshape2 readr ggplot2 smoother stringr tidyr spatstat openssl

ENV LANG=en_US.UTF-8
ENV LC_ADDRESS=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LC_COLLATE=en_US.UTF-8
ENV LC_CTYPE=en_US.UTF-8
ENV LC_IDENTIFICATION=en_US.UTF-8
ENV LC_MEASUREMENT=en_US.UTF-8
ENV LC_MESSAGES=en_US.UTF-8
ENV LC_MONETARY=en_US.UTF-8
ENV LC_NAME=en_US.UTF-8
ENV LC_NUMERIC=en_US.UTF-8
ENV LC_PAPER=en_US.UTF-8
ENV LC_TELEPHONE=en_US.UTF-8
ENV LC_TIME=en_US.UTF-8

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -zxvf bedtools-2.29.1.tar.gz && \
    cd bedtools2 && \
    make && \
    make install

RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -vxjf samtools-1.15.tar.bz2 && \
    cd samtools-1.15 && \
    make && \
    make install

# Copy dependencies into docker image
COPY coverage.R $DESTINATION/coverage.R
COPY checksum.R $DESTINATION/checksum.R
COPY align.R $DESTINATION/align.R
COPY utils.R $DESTINATION/utils.R

RUN ln -s $DESTINATION/coverage.R $DESTINATION/coverage && \
    ln -s $DESTINATION/align.R $DESTINATION/align && \
    ln -s $DESTINATION/checksum.R $DESTINATION/checksum && \
    ln -s $DESTINATION/utils.R /bin/ln_utils.R

RUN sed -i 's/utils.R/\/bin\/ln_utils.R/g' $DESTINATION/align.R && \
    sed -i 's/utils.R/\/bin\/ln_utils.R/g' $DESTINATION/coverage.R && \
    sed -i 's/utils.R/\/bin\/ln_utils.R/g' $DESTINATION/checksum.R

RUN chmod 755 $DESTINATION/coverage $DESTINATION/align $DESTINATION/checksum

WORKDIR /mount
