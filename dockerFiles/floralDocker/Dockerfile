#From public.ecr.aws/lts/ubuntu:20.04_stable
FROM rocker/shiny:4.0.5

MAINTAINER Priyansh Srivastava <psrivastava@bioba.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN echo "Install begin"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
     lzma-dev \
     liblzma-dev \
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
RUN echo "Installed R2"

WORKDIR /

# Installing requiremnts
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shiny", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("tidyverse", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shinythemes", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("shinycssloaders", repos = "https://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("reshape2", repos = "https://cran.us.r-project.org")'

RUN Rscript -e 'BiocManager::install(c("flowCore","ggcyto","flowStats","viridis","flowViz"), update = TRUE, ask = FALSE, Ncpus = 4)'

RUN Rscript -e 'suppressPackageStartupMessages(require(tidyverse)); suppressPackageStartupMessages(require(shiny)); suppressPackageStartupMessages(require(shinythemes));'
RUN Rscript -e 'suppressPackageStartupMessages(require(shinycssloaders)); suppressPackageStartupMessages(require(reshape2)); suppressPackageStartupMessages(require(flowCore));'
RUN Rscript -e 'suppressPackageStartupMessages(require(ggcyto)); suppressPackageStartupMessages(require(flowStats)); suppressPackageStartupMessages(require(viridis));'
RUN Rscript -e 'suppressPackageStartupMessages(require(flowViz))'
RUN Rscript -e 'install.packages("DT", repos = "https://cran.us.r-project.org")'

RUN echo "Installed R1"

RUN git clone https://github.com/spriyansh/esr.git

RUN chmod 777 esr/*

RUN cp /esr/floral/R/app/* /srv/shiny-server/

USER shiny

EXPOSE 3838

CMD R -e "options('shiny.port'=3838,shiny.host='0.0.0.0');experimentdocker::run_app()"

CMD ["/usr/bin/shiny-server"]

RUN echo "--Complete--"

