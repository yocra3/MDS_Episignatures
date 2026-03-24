#' ---------------------------
#'
#' Purpose of script:
#'
#'  Preprocess  GSE129828
#' 
#' ---------------------------


## Load libraries
library(GEOquery)
library(minfi)
library(methylKit)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


## Download GEO data
GSE129828 <- getGEO("GSE129828", GSEMatrix = TRUE)

## Read methylation files
meth_files <- list.files("data/GSE129828/", pattern = "*.txt.gz", full.names = TRUE)

sample_names <- str_extract(basename(meth_files), "^GSM\\d+")

GSE129828_methylkit <- methRead(as.list(meth_files),
              sample.id = as.list(sample_names),
           assembly="hg19",
           context="CpG",
           treatment = rep(1, length(meth_files)),
           mincov = 10
           )

GSE129828_methylkit_filt <- filterByCoverage(GSE129828_methylkit, 
                            lo.count=10, lo.perc=NULL,
                            hi.count=NULL, hi.perc=99.9)


GSE129828_methylkit_united <- methylKit::unite(GSE129828_methylkit_filt)

save(GSE129828_methylkit_united, file = "results/preprocessing/GSE129828/GSE129828_methylkit_united.Rdata")
GSE129828_betas <- percMethylation(GSE129828_methylkit_united)/100

GSE129828_GR <- as(GSE129828_methylkit_united,"GRanges")
mcols(GSE129828_GR) <- NULL

## Annotate CpG sites
data(Locations)

Locations_GR <- makeGRangesFromDataFrame(Locations, seqnames.field = "chr", start.field = "pos", end.field = "pos")
overlaps <- findOverlaps(GSE129828_GR, Locations_GR + 100)


## Read masking CpG data
mask <- read_table("data/HM450.hg38.mask.tsv.gz")

mask_cpgs <- mask %>% filter(MASK_general == TRUE) %>% pull(probeID)

## Process phenotype data
pheno_data_mod <- pheno_data %>%
   mutate(BIRTH_DATE = as.Date(`Date of birth (dd.mm.yyyy)`, format = "%Y-%m-%d"),
   DIAGNOSIS_DATE = as.Date(`Diagnosis date (dd.mm.yyyy)`, format = "%Y-%m-%d"),
   WHO = `WHO (1=CMML-1; 2=CMML-2)`,
   FAB = `FAB\r\n(1=MD; 2=MP)`,
   HB = `Hb \r\n(g/dL)`,
   PLT = `Platelets\r\n(x10^9)`,
    WBC = `WBC\r\n(x10^9)`,
    ANC = `ANC\r\n(x10^9)`,
    BM_BLASTS = `BM Blasts\r\n  (%)`,
    Transfusion = `Tranfusion dependency \r\n(0=no; 1=yes; 9=unknown)`,
    Splenomegaly = `Splenomegaly (0=no; 1=yes; 9=unknown)`,
    CPSS = `CPSS \r\n(0=low risk; 1=intermediate-1; 2=intermediate-2; 3=high risk)`,
    OS_STATUS = `STATUS (0=alive; 1=dead)`,
    LAST_VISIT_DATE = as.Date(`Date of last follow-up/death\r\n(dd.mm.yyyy)`, format = "%Y-%m-%d"),
    AML_PROG = `PROGRESSION TO AML \r\n(1=yes; 2=no)`,
    AMLT_DATE = as.Date(`Progression date (dd.mm.yyyy)`, format = "%Y-%m-%d")
    ) %>%
    select(
      -`Date of birth (dd.mm.yyyy)`,
      -`Diagnosis date (dd.mm.yyyy)`,
      -`WHO (1=CMML-1; 2=CMML-2)`,
      -`FAB\r\n(1=MD; 2=MP)`,
      -`Hb \r\n(g/dL)`,
      -`Platelets\r\n(x10^9)`,
      -`WBC\r\n(x10^9)`,
      -`ANC\r\n(x10^9)`,
      -`BM Blasts\r\n  (%)`,
      -`Tranfusion dependency \r\n(0=no; 1=yes; 9=unknown)`,
      -`Splenomegaly (0=no; 1=yes; 9=unknown)`,
      -`CPSS \r\n(0=low risk; 1=intermediate-1; 2=intermediate-2; 3=high risk)`,
      -`STATUS (0=alive; 1=dead)`,
      -`Date of last follow-up/death\r\n(dd.mm.yyyy)`,
      -`PROGRESSION TO AML \r\n(1=yes; 2=no)`,
      -`Progression date (dd.mm.yyyy)`
    ) %>%
    mutate(Age = as.numeric(difftime(DIAGNOSIS_DATE, BIRTH_DATE, units = "weeks")) / 52.25,
    OS_YEARS = as.numeric(difftime(LAST_VISIT_DATE, DIAGNOSIS_DATE, units = "weeks")) / 52.25, 
    AMLt_STATUS = ifelse(AML_PROG == 1, 1, 0),
    AMLt_YEARS = ifelse(AML_PROG == 1, 
        as.numeric(difftime(AMLT_DATE, DIAGNOSIS_DATE, units = "weeks")) / 52.25, 
        OS_YEARS),
    LFS_STATUS = ifelse(OS_STATUS == 1 | AMLt_STATUS == 1, 1, 0),
    LFS_YEARS = pmin(OS_YEARS, AMLt_YEARS)
    )

## Process Karyotype data
karyotype_data <- pheno_data_mod %>%
    select(sample_title, Karyotype) %>%
    mutate(plus8 = ifelse(str_detect(Karyotype, "\\+8"), 1, 0),
           delY = ifelse(str_detect(Karyotype, "-Y"), 1, 0),
           del7 = ifelse(str_detect(Karyotype, "-7"), 1, 0),
           del20q = ifelse(str_detect(Karyotype, "del\\(20\\)\\(q") | str_detect(Karyotype, "20q-"), 1, 0),
           del12p = ifelse(str_detect(Karyotype, "del\\(12\\)\\(p"), 1, 0),
           plus9 = ifelse(str_detect(Karyotype, "\\+9"), 1, 0),
           del7q = ifelse(str_detect(Karyotype, "del\\(7q\\)"), 1, 0)
           )

## Process mutational data
all_genes <- unique(molecular_data$Gene)
all_genes <- all_genes[all_genes != "Nomutaciones"]
mut_mat <- matrix(0, nrow = nrow(pheno_data), ncol = length(all_genes),
                  dimnames = list(pheno_data$sample_title, all_genes))

mut_mat[subset(pheno_data, is.na(Molecular) | Molecular != "See details of NGS panel in sheet #2")$sample_title, ] <- NA

molecular_data_id <- molecular_data %>%
    mutate(CMML_ID = as.character(`Patient ID`)) %>% 
    left_join(pheno_data_mod %>% select(sample_title, CMML_ID), by = "CMML_ID") %>%
    filter(Gene != "Nomutaciones") 

for (patient in unique(molecular_data_id$sample_title)) {
    patient_gene <- molecular_data_id %>%
        filter(sample_title == patient) %>%
        group_by(Gene) %>%
        summarize(VAF = pmin(sum(`VAF (%)`), 100))
    for (gene in patient_gene$Gene) {
        mut_mat[patient, gene] <- patient_gene %>% filter(Gene == gene) %>% pull(VAF) / 100
    }
}


GSE105420_colData <- pheno_data_mod %>%
    left_join(karyotype_data %>% select(-Karyotype), by = "sample_title") %>%
    left_join(as.data.frame(mut_mat) %>% rownames_to_column("sample_title"), by = "sample_title") %>%
    as.data.frame()
rownames(GSE105420_colData) <- GSE105420_colData$GSM_ID

## Create GenomicRatioSet object
GSE105420_grset <- makeGenomicRatioSetFromMatrix(exprs(GSE105420[[1]]), what = "Beta",
    pData = GSE105420_colData)

save(GSE105420_grset, file = "results/preprocessing/GSE105420/GSE105420_grset_full.Rdata")


GSE105420_grset_mask <- GSE105420_grset[!featureNames(GSE105420_grset) %in% mask_cpgs, ]
save(GSE105420_grset_mask, file = "results/preprocessing/GSE105420/GSE105420_grset_masked.Rdata")

## Make dataset with only case samples
GSE105420_cases <- GSE105420_grset_mask[, !is.na(GSE105420_grset_mask$Karyotype)]
save(GSE105420_cases, file = "results/preprocessing/GSE105420/GSE105420_cases.Rdata")