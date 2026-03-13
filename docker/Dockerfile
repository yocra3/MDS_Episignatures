FROM bioconductor/bioconductor_docker:RELEASE_3_22

## Version 1.0
RUN R -e "BiocManager::install(c('remotes', 'tidyverse', 'pheatmap', 'e1071', 'GEOquery'))"
RUN R -e "BiocManager::install('minfi')"
RUN R -e "remotes::install_github('perishky/meffil')"
