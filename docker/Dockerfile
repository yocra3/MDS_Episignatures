FROM bioconductor/bioconductor_docker:RELEASE_3_22

## Version 1.0
RUN R -e "BiocManager::install(c( \
    'remotes', \
    'tidyverse', \
    'pheatmap', \
    'e1071', \
    'GEOquery', \
    'minfi', \
    'curatedTCGAData', \
    'IlluminaHumanMethylation450kanno.ilmn12.hg19', \
    'IlluminaHumanMethylation450kmanifest', \
    'IlluminaHumanMethylationEPICanno.ilm10b4.hg19', \
    'IlluminaHumanMethylationEPICmanifest', \
    'IlluminaHumanMethylationEPICv2anno.20a1.hg38', \
    'IlluminaHumanMethylationEPICv2manifest' \
)); remotes::install_github('perishky/meffil')"