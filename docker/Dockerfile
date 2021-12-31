FROM rocker/tidyverse:4.1.2

LABEL authors="Viktor Petukhov <viktor.s.petuhov@ya.ru>, Evan Biederstedt <evan.biederstedt@gmail.com>" \
    version.image="1.4.4" \
    version.pagoda2="1.4.4" \
    description="tidyverse image R 4.1.2 to run conos with Rstudio"


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


RUN R -e 'chooseCRANmirror(ind=52); install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "scde", "BiocParallel"))'

RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('testthat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('p2data',dependencies=TRUE, repos='https://kharchenkolab.github.io/drat/', type='source')"

RUN R -e "install.packages('pagoda2',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('conosPanel',dependencies=TRUE, repos='https://kharchenkolab.github.io/drat/', type='source')"

RUN R -e 'devtools::install_github("kharchenkolab/conos")'
