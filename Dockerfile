FROM ubuntu:20.04
LABEL description="TYpeSTeR"
LABEL base_image="ubuntu:latest"
LABEL software="graph genotyper"
LABEL about.home="https://github.com/giacomomutti/TYpeSTeR.git"
LABEL about.license="GPLv3"

ARG DEBIAN_FRONTEND=noninteractive
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

##install TRF

RUN git clone https://github.com/Benson-Genomics-Lab/TRF.git
RUN cd TRF \
    && mkdir build \
    && cd build \
    && ../configure \
    && make 

#no need to env path, because of make install - should be sufficient
ENV PATH /opt/TRF/build/src:$PATH

##install HipSTR

RUN git clone https://github.com/HipSTR-Tool/HipSTR
RUN cd HipSTR \
    && make 

ENV PATH /opt/HipSTR:$PATH


##clone repo with commands required to run the entire genotyping + custom python script
RUN git clone https://github.com/giacomomutti/TYpeSTeR.git
RUN cd TYpeSTeR \
    && chmod 775 typester.sh

ENV PATH /opt/TYpeSTeR:$PATH