#' ---------------------------
#'
#' Purpose of script:
#'
#' Download TCGA data
#' 
#' ---------------------------

## Load libraries
library(curatedTCGAData)
library(MultAssayExperiment)
library(tidyverse)
library(TCGAutils)
library(minfi)

## Install TCGA utils only to check this
data('diseaseCodes', package = "TCGAutils")

curatedTCGAData(
    diseaseCode = "LAML", assays = "*", version = "2.0.1"
)

## Download LAML data
laml_mae <- curatedTCGAData(
    "LAML", c("Methylation_methyl450", "Mutation"), version = "2.0.1", dry.run = FALSE
)

sampleTables(laml_mae)

laml_mae <- intersectColumns(laml_mae)

laml_mut <- laml_mae[[1]]
laml_mut_filt <- laml_mut[mcols(laml_mut)$gene_name %in% c("ASXL1", "TET2", "SRSF2", "RUNX1", "CBL"), ]

laml_methy <- laml_mae[[2]]

extractVAF <- function(gr, gene){
        GR <- subset(gr, Hugo_Symbol == gene)
        if (length(GR) == 0) {
            return(0)
        } else {
            return(sum(as.numeric(GR$TumorVAF_WU))/100)
        }
    }
    
laml_mut_GR <- as(laml_mut_filt, "GRangesList")
laml_methy$ASXL1 <- sapply(laml_mut_GR,extractVAF, gene = "ASXL1")
laml_methy$TET2 <- sapply(laml_mut_GR,extractVAF, gene = "TET2")
laml_methy$SRSF2 <- sapply(laml_mut_GR,extractVAF, gene = "SRSF2")
laml_methy$RUNX1 <- sapply(laml_mut_GR,extractVAF, gene = "RUNX1")
laml_methy$CBL <- sapply(laml_mut_GR,extractVAF, gene = "CBL")

laml_methy_gset <- makeGenomicRatioSetFromMatrix(as.matrix(assay(laml_methy)), what="Beta", 
    pData = colData(laml_methy))

save(laml_methy_gset, file = "results/preprocessing/TCGA_LAML/TCGA_LAML_methylation_mutation.Rdata")