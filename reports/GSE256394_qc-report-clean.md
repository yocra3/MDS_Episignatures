# QC report
- study: Illumina methylation data
- author: Analyst
- date: 23 March, 2026

## Parameters used for QC


```
## $colour.code
## NULL
## 
## $control.categories
## NULL
## 
## $sex.outlier.sd
## [1] 5
## 
## $meth.unmeth.outlier.sd
## [1] 3
## 
## $control.means.outlier.sd
## [1] 5
## 
## $detectionp.samples.threshold
## [1] 0.1
## 
## $beadnum.samples.threshold
## [1] 0.1
## 
## $detectionp.cpgs.threshold
## [1] 0.1
## 
## $beadnum.cpgs.threshold
## [1] 0.1
## 
## $snp.concordance.threshold
## [1] 0.9
## 
## $sample.genotype.concordance.threshold
## [1] 0.9
## 
## $detection.threshold
## [1] 0.01
## 
## $bead.threshold
## [1] 3
## 
## $sex.cutoff
## [1] -2
```
## Number of samples

There are 176 samples analysed.

## Sex mismatches

To separate females and males, we use the difference of total median intensity for Y chromosome probes and X chromosome probes. This will give two distinct clusters of intensities. Females will be clustered on the left and males on the right. 
There are 1 sex detection outliers, and 0 sex detection mismatches.


|sample.name |predicted.sex |declared.sex |   xy.diff|status  |
|:-----------|:-------------|:------------|---------:|:-------|
|GSM8096924  |M             |M            | -1.106535|outlier |

This is a plot of the difference between median 
chromosome Y and chromosome X probe intensities ("XY diff").
Cutoff for sex detection was
XY diff = -2. Mismatched samples are shown in red. The dashed lines represent 5 SD from  the mean xy difference. Samples that fall in this interval are denoted as outliers.



![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)


## Methylated vs unmethylated
To explore the quality of the samples, it is useful to plot the median methylation intensity against the median unmethylation intensity with the option to color outliers by group.
There are 1 outliers from the meth vs unmeth comparison.
Outliers are samples whose predicted median methylated signal is
more than 3 standard deviations
from the expected (regression line).


|sample.name | methylated| unmethylated|   resids| methylated.lm| upper.lm| lower.lm|outliers |
|:-----------|----------:|------------:|--------:|-------------:|--------:|--------:|:--------|
|GSM8097023  |   2940.435|     1905.855| -675.071|      3615.506| 4278.545| 2952.466|TRUE     |

This is a plot of the methylation signals vs unmethylated signals



![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)


## Control probe means

There were 0 outliers detected based on deviations from mean values for control probes. The beachip arrays contain control probe which can be used to evaluate the quality of specific sample processing steps (staining, extension,target removal, hybridization, bisulfate conversion etc.). For each step, a plot has been generated which shows the control means for each sample. Outliers are deviations from the mean. Some of the control probe categories have a very small number of probes. See Page 222 in this doc: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/infinium-hd-methylation-guide-15019519-01.pdf. The most important control probes are the bisulfite1 and bisulfite2 control probes. 



The distribution of sample control means are plotted here:



![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)


## Sample detection p-values

To further explore the quality of each sample the proportion of probes that didn't pass the detection pvalue has been calculated.
There were 0 samples
with a high proportion of undetected probes
(proportion of probes with
detection p-value > 0.01
is > 0.1).



Distribution:



![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)


## Sample bead numbers


To further assess the quality of each sample the proportion of probes that didn't pass the number of beads threshold has been calculated.
There were 0 samples
with a high proportion of probes with low bead number
(proportion of probes with
bead number < 3
is > 0.1).



Distribution:



![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


## CpG detection p-values

To explore the quality of the probes, the proportion of samples that didn't pass the detection pvalue threshold has been calculated.
There were 1318
probes with only background signal in a high proportion of samples
(proportion of samples with detection p-value > 0.01
is > 0.1).
Manhattan plot shows the proportion of samples.



![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

## Low number of beads per CpG

To further explore the quality of the probes, the proportion of samples that didn't pass the number of beads threshold has been calculated.
There were 588 CpGs
with low bead numbers in a high proportion of samples
(proportion of samples with bead number < 3
is > 0.1).
Manhattan plot of proportion of samples.



![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

## Cellular composition estimates




This section omitted.

## SNP probe beta values

The array includes snp probes which can be used to identify sample swaps by comparing these genotypes to genotype calls from a genotype array. First you could check the quality of these snp probes before using them for sample quality.
Distributions of SNP probe beta values are used to determine the quality of the snp probe and should show 3 peaks, one for each genotype probability.


![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

## Genotype concordance




This section omitted.

## R session information


```
## R version 4.5.2 (2025-10-31)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Etc/UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] meffil_1.6.0                preprocessCore_1.72.0      
##  [3] SmartSVA_0.1.3              RSpectra_0.16-2            
##  [5] isva_1.9                    JADE_2.0-4                 
##  [7] qvalue_2.42.0               gdsfmt_1.46.0              
##  [9] statmod_1.5.1               quadprog_1.5-8             
## [11] DNAcopy_1.84.0              fastICA_1.2-7              
## [13] lme4_2.0-1                  Matrix_1.7-4               
## [15] multcomp_1.4-30             TH.data_1.1-5              
## [17] survival_3.8-3              mvtnorm_1.3-6              
## [19] markdown_2.0                gridExtra_2.3              
## [21] Cairo_1.7-0                 knitr_1.51                 
## [23] reshape2_1.4.5              plyr_1.8.9                 
## [25] sva_3.58.0                  BiocParallel_1.44.0        
## [27] genefilter_1.92.0           mgcv_1.9-4                 
## [29] nlme_3.1-168                limma_3.66.0               
## [31] sandwich_3.1-1              lmtest_0.9-40              
## [33] zoo_1.8-15                  MASS_7.3-65                
## [35] illuminaio_0.52.0           lubridate_1.9.5            
## [37] forcats_1.0.1               stringr_1.6.0              
## [39] dplyr_1.2.0                 purrr_1.2.0                
## [41] readr_2.2.0                 tidyr_1.3.2                
## [43] tibble_3.3.0                ggplot2_4.0.2              
## [45] tidyverse_2.0.0             readxl_1.4.5               
## [47] minfi_1.56.0                bumphunter_1.52.0          
## [49] locfit_1.5-9.12             iterators_1.0.14           
## [51] foreach_1.5.2               Biostrings_2.78.0          
## [53] XVector_0.50.0              SummarizedExperiment_1.40.0
## [55] MatrixGenerics_1.22.0       matrixStats_1.5.0          
## [57] GenomicRanges_1.62.1        Seqinfo_1.0.0              
## [59] IRanges_2.44.0              S4Vectors_0.48.0           
## [61] GEOquery_2.78.0             Biobase_2.70.0             
## [63] BiocGenerics_0.56.0         generics_0.1.4             
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.5.2             BiocIO_1.20.0            
##   [3] bitops_1.0-9              R.oo_1.27.1              
##   [5] cellranger_1.1.0          XML_3.99-0.22            
##   [7] httr2_1.2.1               lifecycle_1.0.5          
##   [9] Rdpack_2.6.6              edgeR_4.8.2              
##  [11] lattice_0.22-7            base64_2.0.2             
##  [13] scrime_1.3.7              magrittr_2.0.4           
##  [15] yaml_2.3.12               otel_0.2.0               
##  [17] doRNG_1.8.6.3             askpass_1.2.1            
##  [19] DBI_1.3.0                 minqa_1.2.8              
##  [21] RColorBrewer_1.1-3        abind_1.4-8              
##  [23] R.utils_2.13.0            RCurl_1.98-1.17          
##  [25] rappdirs_0.3.3            rentrez_1.2.4            
##  [27] annotate_1.88.0           commonmark_2.0.0         
##  [29] DelayedMatrixStats_1.32.0 codetools_0.2-20         
##  [31] DelayedArray_0.36.0       xml2_1.5.1               
##  [33] tidyselect_1.2.1          farver_2.1.2             
##  [35] beanplot_1.3.1            GenomicAlignments_1.46.0 
##  [37] jsonlite_2.0.0            multtest_2.66.0          
##  [39] tools_4.5.2               Rcpp_1.1.1               
##  [41] glue_1.8.0                SparseArray_1.10.9       
##  [43] xfun_0.56                 HDF5Array_1.38.0         
##  [45] withr_3.0.2               fastmap_1.2.0            
##  [47] boot_1.3-32               rhdf5filters_1.22.0      
##  [49] openssl_2.3.5             litedown_0.9             
##  [51] digest_0.6.39             mime_0.13                
##  [53] timechange_0.4.0          R6_2.6.1                 
##  [55] RSQLite_2.4.6             R.methodsS3_1.8.2        
##  [57] cigarillo_1.0.0           h5mread_1.2.1            
##  [59] data.table_1.18.2.1       rtracklayer_1.70.1       
##  [61] httr_1.4.8                S4Arrays_1.10.1          
##  [63] pkgconfig_2.0.3           gtable_0.3.6             
##  [65] blob_1.3.0                S7_0.2.1                 
##  [67] siggenes_1.84.0           clue_0.3-67              
##  [69] scales_1.4.0              png_0.1-9                
##  [71] reformulas_0.4.4          tzdb_0.5.0               
##  [73] rjson_0.2.23              curl_7.0.0               
##  [75] nloptr_2.2.1              cachem_1.1.0             
##  [77] rhdf5_2.54.1              AnnotationDbi_1.72.0     
##  [79] restfulr_0.0.16           pillar_1.11.1            
##  [81] grid_4.5.2                reshape_0.8.10           
##  [83] vctrs_0.7.1               cluster_2.1.8.1          
##  [85] xtable_1.8-8              evaluate_1.0.5           
##  [87] GenomicFeatures_1.62.0    cli_3.6.5                
##  [89] compiler_4.5.2            Rsamtools_2.26.0         
##  [91] rlang_1.1.7               crayon_1.5.3             
##  [93] rngtools_1.5.2            labeling_0.4.3           
##  [95] nor1mix_1.3-3             mclust_6.1.2             
##  [97] stringi_1.8.7             hms_1.1.4                
##  [99] sparseMatrixStats_1.22.0  bit64_4.6.0-1            
## [101] Rhdf5lib_1.32.0           KEGGREST_1.50.0          
## [103] rbibutils_2.4.1           memoise_2.0.1            
## [105] bit_4.6.0
```
