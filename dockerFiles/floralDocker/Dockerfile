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

# Installing R
ENV R_VERSION=R-4.2.0
WORKDIR /usr/local
RUN wget https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz && \
        tar xvf $R_VERSION.tar.gz && \
        rm $R_VERSION.tar.gz && \
        cd $R_VERSION && \
./configure --with-readline=no --with-x=no && make -j 4 && make install
ENV PATH=$PATH:~/.local/bin/ 
RUN R --version

RUN echo "Installed R"

WORKDIR /

# Installing requiremnts
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shiny", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("tidyverse", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shinythemes", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shinycssloaders", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("reshape2", repos = "https://cran.us.r-project.org")'


RUN mkdir -p james

WORKDIR /james

RUN echo "--Complete--"

