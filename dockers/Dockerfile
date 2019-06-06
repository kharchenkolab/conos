FROM rocker/tidyverse:3.5.1

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libcurl4-openssl-dev \
  libssl-dev \
  libgsl0-dev \
  libeigen3-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev \
  libxt-dev \
  libgtk2.0-dev \
  libcairo2-dev \
  xvfb  \
  xauth \
  xfonts-base \
  libz-dev \
  libhdf5-dev

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("RcppEigen", "urltools", "Rtsne"))' && \
  R -e 'devtools::install_url("https://cran.r-project.org/src/contrib/Archive/fpc/fpc_2.1-11.tar.gz")' && \
  R -e 'devtools::install_github("satijalab/seurat", ref="v2.3.3")' && \
  R -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("GO.db", "org.Hs.eg.db","org.Mm.eg.db", "pcaMethods","DESeq2"), suppressUpdates=TRUE);' && \
  R -e 'devtools::install_github("jlmelville/uwot")' && \
  R -e 'devtools::install_github("hms-dbmi/pagoda2")' && \
  R -e 'devtools::install_github("hms-dbmi/conos", ref="dev")'
