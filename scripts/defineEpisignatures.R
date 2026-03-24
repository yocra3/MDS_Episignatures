#' ---------------------------
#'
#' Purpose of script:
#'
#'  Define Episignatures in GSE105420
#' 
#' ---------------------------

## Load libraries
library(minfi)
library(tidyverse)
library(e1071)
library(limma)
library(cowplot)

## Load preprocessed data
load("results/preprocessing/GSE105420/GSE105420_cases.Rdata")
load("results/preprocessing/GSE256394/GSE256394.cases.Rdata")
load("results/preprocessing/TCGA_LAML/TCGA_LAML_methylation_mutation.Rdata")



## ============================================================
## MODULARIZED FUNCTIONS FOR EPISIGNATURE ANALYSIS
## ============================================================

#' Run limma analysis for a given variable
#' @param grset GenomicRatioSet object
#' @param formula Formula for the model
#' @param coef Coefficient to test
#' @return Data frame with limma results
run_limma_analysis <- function(grset, formula, coef) {
  model <- model.matrix(formula, data = colData(grset))
  fit <- lmFit(getBeta(grset), model) %>% eBayes()
  topTable(fit, coef = coef, number = Inf, adjust.method = "BH")
}

#' Select CpGs based on filtering criteria
#' @param limma_table Output from run_limma_analysis
#' @param p_threshold P-value threshold (default: 1e-4)
#' @param fc_threshold Absolute logFC threshold (optional)
#' @return Vector of selected CpG probe IDs
select_cpgs <- function(limma_table, p_threshold = 1e-4, fc_threshold = NULL) {
  if (!is.null(fc_threshold)) {
    rownames(subset(limma_table, P.Value < p_threshold & abs(logFC) > fc_threshold))
  } else {
    rownames(subset(limma_table, P.Value < p_threshold))
  }
}

#' Compute PCA and prepare data for plotting
#' @param beta_matrix Beta matrix (CpGs x samples)
#' @param phenotype Vector of phenotype values
#' @param phenotype_name Name of the phenotype variable
#' @return Data frame with PC1, PC2, and phenotype
compute_pca_data <- function(beta_matrix, phenotype, phenotype_name) {
  pc <- prcomp(t(beta_matrix))
  pc_mat <- pc$x
  tib <- tibble(
    PC1 = pc_mat[, 1],
    PC2 = pc_mat[, 2],
    !!phenotype_name := phenotype
  )
  list(pca = pc, pca_data = tib)
}


#' Compute PCA and prepare data for plotting
#' @param beta_matrix Beta matrix (CpGs x samples)
#' @param phenotype Vector of phenotype values
#' @param phenotype_name Name of the phenotype variable
#' @return Data frame with PC1, PC2, and phenotype
predict_pca_data <- function(pc_model, beta_matrix, phenotype, phenotype_name) {
  pc_mat <- predict(pc_model, t(beta_matrix))
  tib <- tibble(
    PC1 = pc_mat[, 1],
    PC2 = pc_mat[, 2],
    !!phenotype_name := phenotype
  )
  tib
}

#' Create PCA plot
#' @param pca_data Output from compute_pca_data
#' @param phenotype_name Name of the phenotype variable
#' @param is_continuous Whether phenotype is continuous (default: TRUE)
#' @return ggplot object
plot_pca <- function(pca_data, phenotype_name, is_continuous = TRUE) {
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = !!sym(phenotype_name))) +
    geom_point(size = 3) +
    theme_bw()
  
  if (is_continuous) {
    p <- p + scale_color_viridis_c()
  }
  
  return(p)
}

#' Train SVM and make predictions
#' @param beta_train Beta matrix for training (CpGs x samples)
#' @param phenotype_train Training phenotype values
#' @param beta_test Beta matrix for testing (CpGs x samples)
#' @return List with svm model and predictions
train_predict_svm <- function(beta_train, phenotype_train, beta_test) {
  model <- svm(t(beta_train), phenotype_train)
  predictions <- predict(model, t(beta_test))
  
  list(model = model, predictions = predictions)
}

#' Create SVM prediction plot
#' @param predictions_rep Predictions on replication set
#' @param phenotype_rep Actual phenotype values (replication)
#' @param predictions_disc Predictions on discovery set
#' @param phenotype_disc Actual phenotype values (discovery)
#' @return ggplot object
plot_svm_predictions <- function(predictions_rep, phenotype_rep, predictions_disc, phenotype_disc) {
  pred_df <- tibble(
    Predicted = c(predictions_rep, predictions_disc),
    Actual = c(phenotype_rep, phenotype_disc),
    Set = c(rep("Replication", length(predictions_rep)), rep("Discovery", length(predictions_disc)))
  )
  
  ggplot(pred_df, aes(x = Actual, y = Predicted)) +
    geom_point(size = 3) +
    facet_wrap(~Set) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
}


## ============================================================
## SIGNATURE FOR BM BLASTS
## ============================================================

# Define discovery set (70% random split)
set.seed(27)
disc_indices <- sample(seq_len(ncol(GSE105420_cases)), size = floor(0.7 * ncol(GSE105420_cases)), replace = FALSE)
rep_indices <- setdiff(seq_len(ncol(GSE105420_cases)), disc_indices)

bm_blasts_disc <- GSE105420_cases[, disc_indices]
bm_blasts_rep <- GSE105420_cases[, rep_indices]

# Run limma analysis
bm_table <- run_limma_analysis(bm_blasts_disc, ~ BM_BLASTS + WBC + HB + PLT + Age + gender, "BM_BLASTS")

# Select CpGs
top_bm_cpgs <- select_cpgs(bm_table, p_threshold = 1e-3, fc_threshold = 0.01)
# Extract beta values
bm_epi_disc <- getBeta(bm_blasts_disc)[top_bm_cpgs, ]
bm_epi_rep <- getBeta(bm_blasts_rep)[top_bm_cpgs, ]

# Get phenotype values
bm_pheno_disc <- bm_blasts_disc$BM_BLASTS
bm_pheno_rep <- bm_blasts_rep$BM_BLASTS

# Compute PCA
bm_pca_disc <- compute_pca_data(bm_epi_disc, bm_pheno_disc, "BM_BLASTS")
bm_pca_rep <- predict_pca_data(bm_pca_disc$pca, bm_epi_rep, bm_pheno_rep, "BM_BLASTS")

## Plot PCA
bm_plot_disc <- plot_pca(bm_pca_disc$pca_data, "BM_BLASTS", is_continuous = TRUE) + ggtitle("Discovery")
bm_plot_rep <- plot_pca(bm_pca_rep, "BM_BLASTS", is_continuous = TRUE) + ggtitle("Replication")

# Save PCA plot
png("figures/Episignatures/BM_BLASTS.png", width = 600, height = 400)
plot_grid(bm_plot_disc, bm_plot_rep, ncol = 2)
dev.off()

# Train SVM and make predictions
bm_svm_result <- train_predict_svm(bm_epi_disc, bm_pheno_disc, bm_epi_rep)
bm_svm_model <- bm_svm_result$model
bm_pred_rep <- bm_svm_result$predictions
bm_pred_disc <- predict(bm_svm_model, t(bm_epi_disc))

# Create SVM prediction plot
bm_plot_pred <- plot_svm_predictions(bm_pred_rep, bm_pheno_rep, bm_pred_disc, bm_pheno_disc)

# Save SVM prediction plot
png("figures/Episignatures/BM_BLASTS_SVM_prediction.png", width = 1200, height = 600)
print(bm_plot_pred)
dev.off()

# Compute metrics
bm_cor_rep <- cor(bm_pred_rep, bm_pheno_rep)
cat("BM_BLASTS correlation (replication):", bm_cor_rep, "\n")


## ============================================================
## SIGNATURE FOR ASXL1
## ============================================================

laml_methy_gset_filt <- laml_methy_gset[rowSums(is.na(getBeta(laml_methy_gset))) == 0, ]
com_cpgs2 <- intersect(rownames(GSE105420_cases), rownames(laml_methy_gset_filt))

# Define discovery set (samples with ASXL1.y data)
asxl1_disc_indices <- which(!is.na(GSE105420_cases$ASXL1.y))
asxl1_rep_indices <- which(is.na(GSE105420_cases$ASXL1.y) & GSE105420_cases$ASXL1.x %in% c(0, 1))



asxl1_disc <- GSE105420_cases[com_cpgs2, asxl1_disc_indices]
asxl1_rep <- GSE105420_cases[com_cpgs2, asxl1_rep_indices]
asxl1_rep2 <- laml_methy_gset_filt[com_cpgs2, ]
# Run limma analysis
asxl1_table <- run_limma_analysis(asxl1_disc, ~ ASXL1.y + Age + gender + BM_BLASTS + HB + WBC + PLT, "ASXL1.y")

# Select CpGs
top_asxl1_cpgs <- select_cpgs(asxl1_table, p_threshold = 1e-4, fc_threshold = NULL)

# Extract beta values
asxl1_epi_disc <- getBeta(asxl1_disc)[top_asxl1_cpgs, ]
asxl1_epi_rep <- getBeta(asxl1_rep)[top_asxl1_cpgs, ]
asxl1_epi_rep2 <- getBeta(asxl1_rep2)[top_asxl1_cpgs, ]


# Get phenotype values
asxl1_pheno_disc <- asxl1_disc$ASXL1.y
asxl1_pheno_rep <- asxl1_rep$ASXL1.x
asxl1_pheno_rep2 <- asxl1_rep2$ASXL1

# Compute PCA
asxl1_pca_disc <- compute_pca_data(asxl1_epi_disc, asxl1_pheno_disc, "ASXL1")
asxl1_pca_rep <- predict_pca_data(asxl1_pca_disc$pca, asxl1_epi_rep, asxl1_pheno_rep, "ASXL1")
asxl1_pca_rep2 <- predict_pca_data(asxl1_pca_disc$pca, asxl1_epi_rep2, asxl1_pheno_rep2, "ASXL1")


# Create PCA plots
asxl1_plot_disc <- plot_pca(asxl1_pca_disc$pca_data, "ASXL1", is_continuous = TRUE)
asxl1_plot_rep <- plot_pca(asxl1_pca_rep, "ASXL1", is_continuous = FALSE)
asxl1_plot_rep2 <- plot_pca(asxl1_pca_rep2, "ASXL1", is_continuous = TRUE)

# Save PCA plot
png("figures/Episignatures/ASXL1.png", width = 1700, height = 600)
plot_grid(asxl1_plot_disc + ggtitle("Discovery"), 
    asxl1_plot_rep + ggtitle("Replication"), 
    asxl1_plot_rep2 + ggtitle("Replication (Blood)"), ncol = 3)
dev.off()

# Train SVM and make predictions
asxl1_svm_result <- train_predict_svm(asxl1_epi_disc, asxl1_pheno_disc, asxl1_epi_rep)
asxl1_svm_model <- asxl1_svm_result$model
asxl1_pred_rep <- asxl1_svm_result$predictions
asxl1_pred_disc <- predict(asxl1_svm_model, t(asxl1_epi_disc))
asxl1_pred_rep2 <- predict(asxl1_svm_model, t(asxl1_epi_rep2))

# Create SVM prediction plot
asxl1_plot_pred <- tibble(
    Predicted = c(asxl1_pred_rep2, asxl1_pred_rep, asxl1_pred_disc),
    Actual = as.numeric(c(asxl1_pheno_rep2, asxl1_pheno_rep, asxl1_pheno_disc)),
    Set = c(rep("Replication (Blood)", length(asxl1_pred_rep2)), 
    rep("Replication", length(asxl1_pred_rep)), 
    rep("Discovery", length(asxl1_pred_disc)))
  )  %>%
  mutate(Set = factor(Set, levels = c("Discovery", "Replication", "Replication (Blood)")))

# Save SVM prediction plot
png("figures/Episignatures/ASXL1_SVM_prediction.png", width = 1200, height = 600)
ggplot(asxl1_plot_pred, aes(x = Actual, y = Predicted)) +
    geom_point(size = 3) +
    facet_wrap(~Set) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
dev.off()

# Confusion matrix for ASXL1
table(asxl1_pred_rep > 0.2, asxl1_pheno_rep)
#        asxl1_pheno_rep
#         0 1
#   FALSE 4 1
#   TRUE  0 6


## ============================================================
## SIGNATURE FOR SRSF2
## ============================================================

com_cpgs <- intersect(rownames(GSE105420_cases), rownames(gset))

# Define discovery set (samples with SRSF2.y data)
srsf2_disc <- GSE105420_cases[com_cpgs, !is.na(GSE105420_cases$SRSF2.y)]
srsf2_rep1 <- GSE105420_cases[com_cpgs, is.na(GSE105420_cases$SRSF2.y)]
srsf2_rep2 <- GSE256394_cases[com_cpgs, ]

# Run limma analysis
srsf2_table <- run_limma_analysis(srsf2_disc, ~ SRSF2.y + Age + gender + BM_BLASTS + HB + WBC + PLT, "SRSF2.y")

# Select CpGs
top_srsf2_cpgs <- select_cpgs(srsf2_table, p_threshold = 1e-3, fc_threshold = 0.5)
#top_srsf2_cpgs <- rownames(srsf2_table)[1:1000]


# Extract beta values
srsf2_epi_disc <- getBeta(srsf2_disc)[top_srsf2_cpgs, ]
srsf2_epi_rep1 <- getBeta(srsf2_rep1)[top_srsf2_cpgs, ]
srsf2_epi_rep2 <- getBeta(srsf2_rep2)[top_srsf2_cpgs, ]

# Get phenotype values
srsf2_pheno_disc <- srsf2_disc$SRSF2.y
srsf2_pheno_rep1 <- srsf2_rep1$SRSF2.x
srsf2_pheno_rep2 <- srsf2_rep2$SRSF2

# Compute PCA
srsf2_pca_disc <- compute_pca_data(srsf2_epi_disc, srsf2_pheno_disc, "SRSF2")
srsf2_pca_rep1 <- predict_pca_data(srsf2_pca_disc$pca, srsf2_epi_rep1, srsf2_pheno_rep1, "SRSF2")
srsf2_pca_rep2 <- predict_pca_data(srsf2_pca_disc$pca, srsf2_epi_rep2, srsf2_pheno_rep2, "SRSF2")
srsf2_pca_rep2$SRSF2 <- as.factor(srsf2_pca_rep2$SRSF2)

# Create PCA plots
srsf2_plot_disc <- plot_pca(srsf2_pca_disc$pca_data, "SRSF2", is_continuous = TRUE)
srsf2_plot_rep1 <- plot_pca(srsf2_pca_rep1, "SRSF2", is_continuous = FALSE)
srsf2_plot_rep2 <- plot_pca(srsf2_pca_rep2, "SRSF2", is_continuous = FALSE)

# Save PCA plot
png("figures/Episignatures/SRSF2.png", width = 1600, height = 600)
plot_grid(srsf2_plot_disc + ggtitle("Discovery"), 
    srsf2_plot_rep1 + ggtitle("Replication (Same dataset)"), 
    srsf2_plot_rep2 + ggtitle("Replication (External dataset)"), ncol = 3)
dev.off()

# Train SVM and make predictions
srsf2_svm_result <- train_predict_svm(srsf2_epi_disc, srsf2_pheno_disc, srsf2_epi_rep1)
srsf2_svm_model <- srsf2_svm_result$model
srsf2_pred_rep1 <- srsf2_svm_result$predictions
srsf2_pred_disc <- predict(srsf2_svm_model, t(srsf2_epi_disc))
srsf2_pred_rep2 <- predict(srsf2_svm_model, t(srsf2_epi_rep2))

# Create SVM prediction plot
srsf2_plot_pred <- plot_svm_predictions(srsf2_pred_rep1, as.numeric(srsf2_pheno_rep1), srsf2_pred_disc, srsf2_pheno_disc)
srsf2_pred_df <- tibble(
    Predicted = c(srsf2_pred_rep2, srsf2_pred_rep1, srsf2_pred_disc),
    Actual = as.numeric(c(srsf2_pheno_rep2, srsf2_pheno_rep1, srsf2_pheno_disc)),
    Set = c(rep("Replication (External dataset)", length(srsf2_pred_rep2)), 
    rep("Replication (Same dataset)", length(srsf2_pred_rep1)), 
    rep("Discovery", length(srsf2_pred_disc)))
  ) %>%
  mutate(Set = factor(Set, levels = c("Discovery", "Replication (Same dataset)", "Replication (External dataset)")))



# Save SVM prediction plot
png("figures/Episignatures/SRSF2_SVM_prediction.png", width = 1700, height = 600)
 ggplot(srsf2_pred_df, aes(x = Actual, y = Predicted)) +
    geom_point(size = 3) +
    facet_wrap(~Set) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
dev.off()

# Confusion matrix for SRSF2
table(srsf2_pred_rep1 > 0.2, srsf2_pheno_rep1)

#         0 1 UNK
#   FALSE 7 1   5
#   TRUE  0 3   2

table(srsf2_pred_rep2 > 0.2, srsf2_pheno_rep2)
#           0   1
#   FALSE 131  19
#   TRUE    5  12

## ============================================================
## SIGNATURE FOR U2AF1
## ============================================================


# Define discovery set (samples with U2AF1 data)
u2af1_disc <- GSE256394_cases[com_cpgs, ]
u2af1_rep <- GSE105420_cases[com_cpgs, !is.na(GSE105420_cases$U2AF1)]

# Run limma analysis
u2af1_table <- run_limma_analysis(u2af1_disc, ~ U2AF1 + Sex + Progression + mdsdiag, "U2AF1")

# Select CpGs
top_u2af1_cpgs <- select_cpgs(u2af1_table, p_threshold = 1e-3, fc_threshold = 0.1)

# Extract beta values
u2af1_epi_disc <- getBeta(u2af1_disc)[top_u2af1_cpgs, ]
u2af1_epi_rep <- getBeta(u2af1_rep)[top_u2af1_cpgs, ]

# Get phenotype values
u2af1_pheno_disc <- u2af1_disc$U2AF1
u2af1_pheno_rep <- u2af1_rep$U2AF1

# Compute PCA
u2af1_pca_disc <- compute_pca_data(u2af1_epi_disc, u2af1_pheno_disc, "U2AF1")
u2af1_pca_rep <- predict_pca_data(u2af1_pca_disc$pca, u2af1_epi_rep, u2af1_pheno_rep, "U2AF1")
u2af1_pca_disc$pca_data$U2AF1 <- as.factor(u2af1_pca_disc$pca_data$U2AF1)
# Create PCA plots
u2af1_plot_disc <- plot_pca(u2af1_pca_disc$pca_data, "U2AF1", is_continuous = FALSE)
u2af1_plot_rep <- plot_pca(u2af1_pca_rep, "U2AF1", is_continuous = TRUE)

# Save PCA plot
png("figures/Episignatures/U2AF1.png", width = 1200, height = 600)
print(plot_grid(u2af1_plot_disc + ggtitle("Discovery"), u2af1_plot_rep + ggtitle("Replication")))
dev.off()


## ============================================================
## SIGNATURE FOR TET2
## ============================================================

set.seed(27)
tet2 <- GSE105420_cases[, !is.na(GSE105420_cases$TET2)]
tet2_disc_indices <- sample(seq_len(ncol(tet2)), size = floor(0.7 * ncol(tet2)), replace = FALSE)
tet2_rep_indices <- setdiff(seq_len(ncol(tet2)), tet2_disc_indices)

tet2_disc <- tet2[com_cpgs2, tet2_disc_indices]
tet2_rep <- tet2[com_cpgs2, tet2_rep_indices]    
tet2_rep2 <- laml_methy_gset_filt[com_cpgs2, ]


# Run limma analysis
tet2_table <- run_limma_analysis(tet2_disc, ~ TET2 + Age + gender + BM_BLASTS + HB + WBC + PLT, "TET2")

# Select CpGs
top_tet2_cpgs <- select_cpgs(tet2_table, p_threshold = 1e-4, fc_threshold = NULL)

# Extract beta values
tet2_epi_disc <- getBeta(tet2_disc)[top_tet2_cpgs, ]
tet2_epi_rep <- getBeta(tet2_rep)[top_tet2_cpgs, ]
tet2_epi_rep2 <- getBeta(tet2_rep2)[top_tet2_cpgs, ]

# Get phenotype values
tet2_pheno_disc <- tet2_disc$TET2
tet2_pheno_rep <- tet2_rep$TET2
tet2_pheno_rep2 <- tet2_rep2$TET2

# Compute PCA
tet2_pca_disc <- compute_pca_data(tet2_epi_disc, tet2_pheno_disc, "TET2")
tet2_pca_rep <- predict_pca_data(tet2_pca_disc$pca, tet2_epi_rep, tet2_pheno_rep, "TET2")
tet2_pca_rep2 <- predict_pca_data(tet2_pca_disc$pca, tet2_epi_rep2, tet2_pheno_rep2, "TET2")

# Create PCA plots
tet2_plot_disc <- plot_pca(tet2_pca_disc$pca_data, "TET2", is_continuous = TRUE)
tet2_plot_rep <- plot_pca(tet2_pca_rep, "TET2", is_continuous = TRUE)
tet2_plot_rep2 <- plot_pca(tet2_pca_rep2, "TET2", is_continuous = TRUE)

# Save PCA plot
png("figures/Episignatures/TET2.png", width = 1700, height = 600)
plot_grid(tet2_plot_disc + ggtitle("Discovery"), 
    tet2_plot_rep + ggtitle("Replication"), 
    tet2_plot_rep2 + ggtitle("Replication (Blood)"),
    ncol = 3)
dev.off()

# Train SVM and make predictions
tet2_svm_result <- train_predict_svm(tet2_epi_disc, tet2_pheno_disc, tet2_epi_rep)
tet2_svm_model <- tet2_svm_result$model
tet2_pred_rep <- tet2_svm_result$predictions
tet2_pred_disc <- predict(tet2_svm_model, t(tet2_epi_disc))
tet2_pred_rep2 <- predict(tet2_svm_model, t(tet2_epi_rep2))

# Create SVM prediction plot
tet2_plot_pred <- plot_svm_predictions(tet2_pred_rep, as.numeric(tet2_pheno_rep), tet2_pred_disc, tet2_pheno_disc)
tet2_pred_df <- tibble(
    Predicted = c(tet2_pred_rep2, tet2_pred_rep, tet2_pred_disc),
    Actual = as.numeric(c(tet2_pheno_rep2, tet2_pheno_rep, tet2_pheno_disc)),
    Set = c(rep("Replication (Blood)", length(tet2_pred_rep2)), 
    rep("Replication", length(tet2_pred_rep)), 
    rep("Discovery", length(tet2_pred_disc)))
  ) %>%
  mutate(Set = factor(Set, levels = c("Discovery", "Replication", "Replication (Blood)")))



# Save SVM prediction plot
png("figures/Episignatures/TET2_SVM_prediction.png", width = 1700, height = 600)
 ggplot(tet2_pred_df, aes(x = Actual, y = Predicted)) +
    geom_point(size = 3) +
    facet_wrap(~Set) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
dev.off()

cor(tet2_pred_rep, tet2_pheno_rep)
cor(tet2_pred_rep2, tet2_pheno_rep2)

table(tet2_pred_rep > 0.2, tet2_pheno_rep > 0)

## ============================================================
## SIGNATURE FOR RUNX1
## ============================================================

# Select data
runx1_disc <- GSE105420_cases[com_cpgs2, !is.na(GSE105420_cases$RUNX1) & GSE105420_cases$RUNX1 < 0.5]
runx1_rep <-  laml_methy_gset_filt[com_cpgs2, ]

# Run limma analysis
runx1_table <- run_limma_analysis(runx1_disc, ~ RUNX1 + Age + gender + BM_BLASTS + HB + WBC + PLT, "RUNX1")

# Select CpGs
top_runx1_cpgs <- select_cpgs(runx1_table, p_threshold = 1e-4, fc_threshold = NULL)

# Extract beta values
runx1_epi_disc <- getBeta(runx1_disc)[top_runx1_cpgs, ]
runx1_epi_rep <- getBeta(runx1_rep)[top_runx1_cpgs, ]

# Get phenotype values
runx1_pheno_disc <- runx1_disc$RUNX1
runx1_pheno_rep <- runx1_rep$RUNX1

# Compute PCA
runx1_pca_disc <- compute_pca_data(runx1_epi_disc, runx1_pheno_disc, "RUNX1")
runx1_pca_rep <- predict_pca_data(runx1_pca_disc$pca, runx1_epi_rep, runx1_pheno_rep, "RUNX1")

# Create PCA plots
runx1_plot_disc <- plot_pca(runx1_pca_disc$pca_data, "RUNX1", is_continuous = TRUE)
runx1_plot_rep <- plot_pca(runx1_pca_rep, "RUNX1", is_continuous = TRUE)

# Save PCA plot
png("figures/Episignatures/RUNX1.png", width = 1200, height = 600)
print(plot_grid(runx1_plot_disc + ggtitle("Discovery"), runx1_plot_rep + ggtitle("Replication (Blood)")))
dev.off()

## Train SVM and make predictions
runx1_svm_result <- train_predict_svm(runx1_epi_disc, runx1_pheno_disc, runx1_epi_rep)
runx1_svm_model <- runx1_svm_result$model
runx1_pred_rep <- runx1_svm_result$predictions
runx1_pred_disc <- predict(runx1_svm_model, t(runx1_epi_disc))

# Create SVM prediction plot
runx1_plot_pred <- plot_svm_predictions(runx1_pred_rep, as.numeric(runx1_pheno_rep), runx1_pred_disc, runx1_pheno_disc)

# Save SVM prediction plot
png("figures/Episignatures/RUNX1_SVM_prediction.png", width = 1200, height = 600)
print(runx1_plot_pred)
dev.off()