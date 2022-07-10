From public.ecr.aws/lts/ubuntu:20.04_stable

MAINTAINER Priyansh Srivastava <psrivastava@bioba.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN echo "Install begin"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libmysqlclient-dev \
    libpng-dev \
    libmysqlclient-dev \
    libncurses-dev \
    libssl-dev \
    python3 \
    python3-pip \
    python3-dev \
    gfortran \
    perl \
    liblzma-dev \
    libpcre2-dev \
    libcurl4-openssl-dev \
    libreadline-gplv2-dev \
    libncursesw5-dev \
    libncurses5-dev \
    libssl-dev \
    libsqlite3-dev \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    libbz2-dev \
    libffi-dev \
    zlib1g-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    wget \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    libcairo2-dev \
    libxt-dev \
    git-all

RUN apt-get update

WORKDIR /

RUN mkdir -p james

WORKDIR /james

# Break dancer installation begins
RUN git clone --recursive https://github.com/genome/breakdancer.git
RUN echo "Cloning Success"
WORKDIR /james/breakdancer
RUN echo "In breakdancer"
RUN mkdir -p build
RUN echo "Created build"
WORKDIR /james/breakdancer/build
RUN echo "In Build build"
RUN cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
RUN echo "Cmake done"
RUN make
RUN echo "make done"
RUN make install
WORKDIR /james
RUN echo "success"
RUN apt install -y cpanminus libgd-perl libgd-dev libgd-graph-perl
RUN cpan install GD 
RUN cpanm Statistics::Descriptive
RUN cpanm GD::Graph::histogram
RUN apt-get install -y samtools

RUN perl breakdancer/perl/bam2cfg.pl breakdancer/test-data/NA19238_chr21_del_inv.bam Shit_it_runs.cfg

