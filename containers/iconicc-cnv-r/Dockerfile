FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/2025-01-15"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

RUN Rscript -e 'install.packages(c( \
      "data.table", "reshape2", "scales", "latex2exp", "ggplot2", "useful", \
      "tidyverse", "pracma", "optparse", "dplyr", "gridExtra", "Hmisc", "remotes" \
    ))'

RUN Rscript -e 'BiocManager::install(c("rhdf5", "rCGH", "DNAcopy", "GenomicRanges"), \
      update = FALSE, ask = FALSE)'

RUN Rscript -e 'remotes::install_github( \
      "mskilab-org/gUtils@095377ab04195aecf6ebb8397072b24864c81599")'
