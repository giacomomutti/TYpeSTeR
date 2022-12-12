FROM ubuntu:20.04
LABEL description="TYpeSTeR"
LABEL base_image="ubuntu:latest"
LABEL software="graph genotyper"
LABEL about.home="https://github.com/giacomomutti/TYpeSTeR.git"
LABEL about.license="GPLv3"

ARG DEBIAN_FRONTEND=noninteractive
ARG R_VERSION=3.5.2

#install basic libraries and python
WORKDIR /opt

RUN apt-get update
RUN apt-get -y install build-essential \
	wget git\
	bzip2 libbz2-dev \
	zlib1g zlib1g-dev \
	liblzma-dev \
	libssl-dev \
	libncurses5-dev \
	libz-dev \
	python3-distutils python3-dev python3-pip \ 
	libjemalloc-dev \
	cmake make g++ \
	libhts-dev \
	libzstd-dev \
	autoconf \
	libatomic-ops-dev \
	pkg-config \
	cargo \
	pigz \
    	&& apt-get -y clean all \
    	&& rm -rf /var/cache

RUN ln -s /usr/bin/python3 /usr/bin/python


# # Install system dependencies for the tidyverse R packages
# RUN apt-get install -y \
#     make \
#     libcurl4-openssl-dev \
#     libssl-dev \
#     pandoc \
#     libxml2-dev \
#     gdebi-core \
#     	&& apt-get -y clean all \
#     	&& rm -rf /var/cache

# # # Install R
# # # download a version of R and build from source
# # RUN wget https://cdn.rstudio.com/r/ubuntu-1604/pkgs/r-${R_VERSION}_1_amd64.deb
# # RUN gdebi r-${R_VERSION}_1_amd64.deb


# Install R packages
#  RUN R -e "install.packages(c('methods', 'jsonlite', 'tseries'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
	&& tar -jxvf samtools-1.16.1.tar.bz2 \
	&& rm samtools-1.16.1.tar.bz2 \
	&& cd samtools-1.16.1 \
	&& ./configure \
	&& make \
	&& make install
#no need to env path, because of make install - should be sufficient

# install bedtools
RUN mkdir -p bedtools \
	&& cd bedtools \
	&& wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
	&& mv bedtools.static.binary bedtools \
	&& chmod a+x bedtools
ENV PATH /opt/bedtools:$PATH

# Install blast+ tool

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz \
	&& tar xvfz ncbi-blast-2.13.0+-x64-linux.tar.gz \
	&& rm ncbi-blast-2.13.0+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-2.13.0+/bin:$PATH


##install TRF
RUN git clone https://github.com/Benson-Genomics-Lab/TRF.git
RUN cd TRF \
    && mkdir build \
    && cd build \
    && ../configure \
    && make 

ENV PATH /opt/TRF/build/src:$PATH

##install HipSTR

RUN git clone https://github.com/HipSTR-Tool/HipSTR
RUN cd HipSTR \
    && make 

ENV PATH /opt/HipSTR:$PATH

##clone repo with commands required to run the entire genotyping + custom python script 
RUN git clone https://github.com/giacomomutti/TYpeSTeR.git
RUN cd TYpeSTeR \
    && chmod 775 typester.sh \
    && chmod 775 scripts/find_homology.sh

ENV PATH /opt/TYpeSTeR:$PATH
