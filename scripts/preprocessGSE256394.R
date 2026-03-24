#' ---------------------------
#'
#' Purpose of script:
#'
#'  Preprocess  GSE256394
#' 
#' wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256394/suppl/GSE256394_RAW.tar
#' ---------------------------


## Load libraries
library(GEOquery)
library(minfi)
library(readxl)
library(tidyverse)
library(meffil)


## Download GEO data
GSE256394 <- getGEO("GSE256394", GSEMatrix = TRUE)

pheno <- pData(GSE256394[[1]])

outPrefix <- "results/preprocessing/GSE256394/GSE256394"
## Preprocess with meffil
options(mc.cores=8)


samplesheet <- meffil.create.samplesheet("data/GSE256394/", recursive=TRUE)
samplesheet$Sex <- ifelse(pheno[samplesheet$Sample_Name, "gender:ch1"] == "Male", "M", "F")

qc.objects <- meffil.qc(samplesheet, verbose=TRUE)
save(qc.objects, file = "results/preprocessing/GSE256394/GSE256394.qc.objects.Robj")


qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5
)
qc.summary <- meffil.qc.summary(
	qc.objects,
	parameters = qc.parameters)

meffil.qc.report(qc.summary, output.file="reports/GSE256394_qc-report.html")

outlier <- qc.summary$bad.samples
table(outlier$issue)
criteria <- c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")
index <- outlier$issue %in% criteria
outlier <- outlier[index,]

qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
save(qc.objects, file = "results/preprocessing/GSE256394/GSE256394.qc.objects.clean.Robj")

qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
meffil.qc.report(qc.summary, output.file = "reports/GSE256394_qc-report-clean.html")


outlier2 <- qc.summary$bad.samples
outlier2 <- subset(outlier2, issue %in% criteria)
qc.objects <- meffil.remove.samples(qc.objects, outlier2$sample.name)
save(qc.objects, file = "results/preprocessing/GSE256394/GSE256394.qc.objects.clean2.Robj")
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
meffil.qc.report(qc.summary, output.file = "reports/GSE256394_qc-report-clean2.html")

outs <- rbind(outlier, outlier2)


## Report filtered samples and probes
write.table(outs, file = "results/preprocessing/GSE256394/GSE256394.removed.samples.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = "results/preprocessing/GSE256394/GSE256394.removed.probes.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")


## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = paste0(outPrefix, ".pc.fit.pdf"), height = 6, width = 6)


pcs <- 4
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)

## Add predicted sex as sample sheet variable
for (i in seq_len(length(norm.objects))){
  norm.objects[[i]]$samplesheet$pred.sex <- norm.objects[[i]]$predicted.sex
}
save(norm.objects, file = paste0(outPrefix, ".norm.obj.pc.Rdata"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name,
                                      verbose = TRUE)
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 40000)

batch_var <- c("Slide", "sentrix_col",  "sentrix_row", "Sex")

norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(8),
  batch.pcs = seq_len(8),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.beta, file = paste0(outPrefix, ".norm.beta.Rdata"))
save(norm.summary, file = paste0(outPrefix, ".norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.normalization.html"))

rownames(samplesheet) <- samplesheet$Sample_Name
samplesheet.final <- samplesheet[colnames(norm.beta), ]

## Add cell counts
pheno_final <- left_join(samplesheet.final, pheno, by = c("Sample_Name" = "geo_accession")) %>%
    mutate(Progression = `progression:ch1`,
    mdsdiag = `mdsdiag:ch1`,
    disease = `disease:ch1`)

genes_of_interest <- c("SF3B1", "SRSF2", "U2AF1")

for (gene in genes_of_interest) {
    gene_status <- str_match(
        pheno_final$`genotype:ch1`,
        paste0("\\b", gene, "\\b\\s*[:=-]?\\s*(WT|MUT)\\b")
    )[, 2]

    pheno_final[[gene]] <- case_when(
        gene_status == "WT" ~ 0L,
        gene_status == "MUT" ~ 1L,
        TRUE ~ NA_integer_
    )
}



## Save genomicratioset
GSE256394_all <- makeGenomicRatioSetFromMatrix(norm.beta, pData = pheno_final,
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")
save(GSE256394_all, file = paste0(outPrefix, ".GenomicRatioSet.Rdata"))

## Create GRSets  ####
mask <- read_table("data/EPIC.hg38.mask.tsv.gz")

### Probes not measuring methylation
GSE256394_cg <- dropMethylationLoci(GSE256394_all)
save(GSE256394_cg, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
GSE256394_filt <- GSE256394_cg[!rownames(GSE256394_cg) %in% subset(mask, MASK_general == TRUE)$probeID, ]
save(GSE256394_filt, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
GSE256394_autosomal <- GSE256394_filt[!seqnames(rowRanges(GSE256394_filt)) %in% c("chrX", "chrY"), ]
save(GSE256394_autosomal, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

GSE256394_cases <- GSE256394_autosomal[, GSE256394_autosomal$disease == "LR-MDS"]
save(GSE256394_cases, file = paste0(outPrefix, ".cases.Rdata"))

# ## Create Initial and final dataset with missings
# dp <- meffil.load.detection.pvalues(qc.objects)
# gset <- ori
# dp.f <- dp[rownames(gset), colnames(gset)]

# beta <- assay(gset)
# beta[dp.f > 2e-16] <- NA
# assay(gset) <- beta
# save(gset, file = paste0(outPrefix, ".allCpGs.withNA.GenomicRatioSet.Rdata"))

# gset <- gset[rownames(final), colnames(final)]
# save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

