#===
#Vignette using DCC and PKC Files of NanoString: https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
#Vignette using CountData of NanoString: https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
#===

#### Define Functions########################
# function to create project folder if not same as R Project folder and io folders in it 
#...project folder could be the R project or within it with the script & io within the project
setupProject <- function(project) {
  rstudio_dir <- rstudioapi::getActiveProject() # Get the Rstudio Directory
  rstudio_base <- basename(rstudio_dir) # Get the base name of the Rstudio Directory
  if (rstudio_base != project) { # If the base name != project, create a folder named as the project inside rstudio_dir
    project_dir <- file.path(rstudio_dir, project)
    if (!file.exists(project_dir)) {
      dir.create(project_dir, recursive = TRUE)
    }
  } else {
    # If the base name is the same as the project, use the rstudio_dir as the project_dir
    project_dir <- rstudio_dir
  }
  # Define input and output directories inside the project directory
  input_dir <- file.path(project_dir, "input")
  output_dir <- file.path(project_dir, "output")
  # Check if input and output directories exist, if not, create them
  if (!file.exists(input_dir)) {
    dir.create(input_dir, recursive = TRUE)
  }
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Set the working directory to the project directory
  setwd(project_dir)
  # Print Confirmation of Correct Folder Structure
  if (basename(getwd()) == project) {
    print("Folder setup correctly")
  } else {
    print("Fix folder structure")
  }
  # Export input_dir and output_dir to the global environment
  assign("input_dir", input_dir, envir = .GlobalEnv)
  assign("output_dir", output_dir, envir = .GlobalEnv)
}
savePDF <- function(figure, fileName, w = 7, h = 10) {
  currentDate <- format(Sys.Date(), "%Y%m%d") #current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".pdf"))
  # Save the figure
  pdf(file = fullFilePath, width = w, height = h, pointsize = 300 / 72)
  print(figure)
  dev.off()
}
savePNG <- function(figure, fileName, w = 900, h = 1300) {
  currentDate <- format(Sys.Date(), "%Y%m%d") # current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".png"))
  # Save the figure as PNG with dimensions in pixels
  png(file = fullFilePath, width = w, height = h, units = "px")
  # Render the plot
  print(figure)
  dev.off()
}

#### Analysis Specific Code ----
#Initiate project
setupProject("NanostringAnalysis")
getwd()
# Project specific override: output folder if any
#output_dir <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/RNASeqAnalysis/output/1.0_QualityControl"

#### Install & Load Packages ----
list.of.packages.cran <- c("DT", "ggalluvial", "ggrepel", "igraph", "tidyverse")
new.packages.cran <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages.cran)>0) install.packages(new.packages.cran)
# Install not-yet-installed Bioconductor packages
list.of.packages.bioc <- c("edgeR", "GSEABase", "limma", "msigdb", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0)if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(new.packages.bioc, update = FALSE)
# Load packages
sapply(c(list.of.packages.cran, list.of.packages.bioc), require, character.only=TRUE)

#rm(list=ls()) #remove all existing lists; DONT

#### Source & Process Input files ----
# sampleAnnoFile from SegmentProperties sheet
# countFile & featureAnnoFile from BioProbeCountMatrix sheet

#sampleAnnoFile ----
#SegmentProperties Worksheet: SegmentDisplayName is default column
sampleAnnoFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                                     sheet="SegmentProperties")

final_sampleAnnoFile <- sampleAnnoFile %>% 
  mutate(SlideName = gsub(" breast", "", 
                          gsub(" 1/2-", "", SlideName)),
         ScanLabel = gsub(" 1/2-", "", ScanLabel),
         SegmentDisplayName = gsub(" 1/2-", "", SegmentDisplayName),
         tissue = gsub(" ", "_", tissue)) %>% 
  mutate_at(6:18, ~ as.logical(.)) %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

# write.table(sampleAnnoFile, "output/sampleAnnoFile.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
# duplicates <- c('D830030K20Rik', 'Gm10406', 'LOC118568634')

#featureAnnoFile----
#BioProbeCountMatrix Worksheet: TargetName is default column
featureAnnoFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                                      sheet = "BioProbeCountMatrix")

featureAnnoFile <- dplyr::select(featureAnnoFile, 1:12) %>% 
  dplyr::select(TargetName, everything())

# Isolate rows with "NegProbe-WTX" in TargetName
featureAnnoFile_NegProbeWTX <- featureAnnoFile %>%
  filter(TargetName == "NegProbe-WTX")

# Process other rows to remove duplicates based on "TargetName" and concatenate unique values for duplicates
  # featureAnnoFile_others <- featureAnnoFile %>%
  #   filter(TargetName != "NegProbe-WTX") %>% #Remaining rows with genes in TargetName
  #   mutate_at(c("GeneID"), as.character) %>% #Convert GeneID num -> character
  #   dplyr::group_by(TargetName) %>% #Concatenate unique values of duplicate TargetNames
  #   summarise(across(
  #     .cols = everything(),
  #     .fns = ~ if (n() > 1) {
  #       unique_values <- unique(trimws(unlist(strsplit(as.character(.), ","))))
  #       paste(unique_values, collapse = ",")
  #     } else {
  #       first(.)
  #     }
  #   ), .groups = "drop") %>%
  #   mutate(GeneID = first(featureAnnoFile$GeneID[match(TargetName, featureAnnoFile$TargetName)]))

#Improved way addresses issues in above code.
#.fns = ~ if (n() > 1) { ... } else { .[1] } to handle cases where there are duplicates or not, 
#ensuring .fns works with character vectors properly.
#Added rowwise() before the mutate step to ensure that mutate works row-wise
#Used indexing [1] instead of first to avoid issues with the first method for numeric vectors.
featureAnnoFile_others <- featureAnnoFile %>%
  filter(TargetName != "NegProbe-WTX") %>% # Remaining rows with genes in TargetName
  mutate(GeneID = as.character(GeneID)) %>% # Convert GeneID num -> character
  group_by(TargetName) %>% # Concatenate unique values of duplicate TargetNames
  summarise(across(
    .cols = everything(),
    .fns = ~ if (n() > 1) {
      unique_values <- unique(trimws(unlist(strsplit(as.character(.), ","))))
      paste(unique_values, collapse = ",")
    } else {
      .[1]
    }
  ), .groups = "drop") %>%
  rowwise() %>%
  mutate(GeneID = featureAnnoFile$GeneID[match(TargetName, featureAnnoFile$TargetName)][1])




# Combine both processed datasets
final_featureAnnoFile <- bind_rows(featureAnnoFile_NegProbeWTX, featureAnnoFile_others) %>% 
  as.data.frame(., row.names = NULL, optional = FALSE, stringsAsFactors = FALSE)

#
# Previous attempts
# featureAnnoFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
#                               sheet="BioProbeCountMatrix") %>% 
#   select(., 1:12) %>% 
#   # filter(TargetName != duplicates) %>%
#   select(TargetName, everything()) %>% 
#   group_by(TargetName) %>%
#   # summarise(across(where(~first(TargetName) != "NegProbe-WTX"), ~ {
#   summarise(across(.cols = !TargetName == "NegProbe-WTX", .fns = ~ {
#     if (n() > 1) {
#       unique_values <- unique(trimws(unlist(strsplit(as.character(.), ","))))
#       paste(unique_values, collapse = ",")
#     } else {
#       .
#     }
#   }), .groups = "drop") %>% 
#   as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#====
# duplicated_values <- unique(featureAnnoFile$TargetName[duplicated(featureAnnoFile$TargetName)])

# column_names <- colnames(featureAnnoFile)
# 
# # Check if "TargetName" is in the column names
# if ("TargetName" %in% column_names) {
#   print("Column 'TargetName' exists in the dataset.")
# } else {
#   print("Column 'TargetName' does not exist in the dataset.")
# }

# featureAnnoFile <- featureAnnoFile %>% #filter(!(TargetName %in% c('D830030K20Rik', 'Gm10406', 'LOC118568634'))) %>% 
#   # filter(TargetName != "NegProbe-WTX") %>%
#   group_by(TargetName) %>% 
#   summarise(across(everything(), ~ paste(unique(.), collapse = "|"))) %>% 
#   select(ProbeName, ProbeDisplayName, everything()) %>% 
#   mutate_at(1, as.integer) %>% 
#   as.data.frame()
# 
# is.data.frame(featureAnnoFile)

#CountFile----
countFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                        sheet="BioProbeCountMatrix") 
countFile <- countFile %>%   dplyr::select(., c(3, 13:187)) %>% 
  # filter(TargetName !=duplicates) %>%
  setNames(c(gsub(" 1/2-", "", colnames(.)))) %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#countFile_NegProbeWTX <- countFile[which(countFile[,"TargetName"]=="NegProbe-WTX"),]
countFile_NegProbeWTX <- countFile %>% filter(TargetName == "NegProbe-WTX")

countFile_others <- countFile %>%
  filter(TargetName != "NegProbe-WTX") %>%
  # mutate_at(c("GeneID"), as.character) %>% 
  group_by(TargetName) %>%
  summarise(across(everything(), ~ sum(., na.rm = TRUE)), .groups = "drop")
# Combine
final_countFile <- bind_rows(countFile_NegProbeWTX, countFile_others) %>% 
  as.data.frame(., row.names = NULL, optional = FALSE, stringsAsFactors = FALSE)

#----
#colnames(countFile) %>% print() #All column names
#colnames(countFile) = gsub(" 1/2-", "", colnames(countFile))
#Also known as CountFile: TargetName is default column

# rownames(countFile) <- NULL
# countFile <- countFile %>% #filter(!(TargetName %in% c('D830030K20Rik', 'Gm10406', 'LOC118568634'))) %>% 

#Remove Duplicates in countFile
# countFile <- countFile %>% 
#   filter(TargetName != "NegProbe-WTX") %>%
#   group_by(TargetName) %>% 
#   # summarise(across(everything(), sum, na.rm = TRUE)) %>% 
#   summarise(across(everything(), ~ paste(unique(.), collapse = "|"))) %>% 
#   mutate(across(-TargetName, as.integer)) %>% 
#   as.data.frame() #%>%


# as.matrix()
# rownames(countFile) <- NULL
# rownames(featureAnnoFile) <- NULL
# rownames(featureAnnoFile) %>% print() 
# is.data.frame(countFile)

#Write countFile to output directory (Optional)
# write.table(countFile, "output/countFile.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
#Write featureAnnoFile to output directory
# write.table(featureAnnoFile, "output/featureAnnoFile.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
#Write sampleAnnoFile to output directory
# write.table(sampleAnnoFile, "output/sampleAnnoFile.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# rownames(countFile
# head(countFile)[,1:5]
# head(sampleAnnoFile)[,1:5]
# head(featureAnnoFile)[,1:5]
# 
# colnames(countFile) %>% print()
# colnames(sampleAnnoFile) %>% print()
# colnames(featureAnnoFile) %>% print()
# rownames(featureAnnoFile) %>% print()

# BiocManager::install("standR")
library(standR)
#install.packages("ggalluvial") #required for the standR package
#install.packages("magick") #required for the standR package
# Create spatialExperiement object

seo <- readGeoMx(countFile = final_countFile,
                  sampleAnnoFile = final_sampleAnnoFile,
                  featureAnnoFile = final_featureAnnoFile)

# seo_df <- readGeoMx(countFile = "output/countFile.txt",
                 # sampleAnnoFile = "output/sampleAnnoFile.txt",
                 # featureAnnoFile = "output/featureAnnoFile.txt",
                 # # colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"), #optional as we are using standard
                 # rmNegProbe = TRUE,
                 # NegProbeName = "NegProbe-WTX")

library(SpatialExperiment)
seo

#Examine the data
assays(seo)
seo$ROILabel
seo$AOINucleiCount
names(seo@metadata)
names(seo_qc@metadata)

assayNames(seo)

##View the count table using the assay function
assay(seo, "counts")[1:5, 1:5]
assay(seo, "logcounts")[1:5, 1:5]
#Sample metadata is in the colData object
colData(seo)[1:5,1:5]
#Gene metadata stored in the rowData of the object
rowData(seo)[1:5,1:5]
#WTA has NegProbe-WTX data. Ensure that there are no duplicate gene names in the TargetName column
metadata(seo)$NegProbes[,1:5]
colData(seo)$QCFlags
colData(seo)$tissue
names(seo@colData)

#QC Steps 
#Sample level QC
library(ggplot2)
library(ggalluvial)
plotSampleInfo(seo, column2plot =c("tissue", "SlideName", "CD45", "Neuron"))

#Gene level QC
seo #Dim 19962x175
seo_qc <- addPerROIQC(seo, 
                      rm_genes =TRUE,
                      sample_fraction = 0.9, #Default
                      min_count = 5) #Default
seo_qc #Gene with low count and expression values in more than threshold (sample_fraction=0.9)
# are removed by applying the function. Dim 19948x175, so removed 14 genes
dim(seo) 
dim(seo_qc) #Genes not meeting the above criteria were removed 19962 > 19948 
names(seo_qc@colData)
length(seo_qc@metadata$genes_rm_rawCount)
names(seo_qc@metadata)
metadata(seo_qc) |> names() #Same as above
metadata(seo) |>names()

#addPerROIQC added columns to colData :lib_size, countOfLowEprGene, percentOfLowEprGene
# and also added columns to metadata: lcpm_threshold, genes_rm_rawCount, genes_rm_logCPM  


# plotGeneQC(seo_qc, ordannots = "regions", col = regions, point_size = 2)
plotGeneQC(seo_qc)
plotGeneQC(seo_qc, top_n=12, ordannots = "tissue", col = tissue, point_size = 2)

#colData(seo)$regions
#colnames(seo) %>% print() 

# data("seo_qc")
# seo_subset <-  addPerROIQC(seo_qc) #Wrong; already done this step to get seo_qc
# plotGeneQC(seo_subset)


#ROI level QC
plotROIQC(seo_qc, x_threshold = 150, color = SlideName)

colData(seo_qc)$AOINucleiCount 
#same as
seo_qc@colData$AOINucleiCount

#AOINuclei count of 150 looks like a good threshold from the figure
qc <- colData(seo_qc)$AOINucleiCount > 150
table(qc) # 3 Values Below threshold
dim(seo_qc) # Dim 19948x175
seo_qc_roi <- seo_qc[, qc]
dim(seo_qc_roi) # We removed 3 ROI (samples/columns)from dataset. Dim 19948x172
# Comparing the Library Size with ROI Area size

plotROIQC(seo_qc_roi, 
          x_threshold = 20000, 
          x_axis = "AOISurfaceArea", 
          x_lab = "AreaSize", 
          y_axis = "lib_size", 
          y_lab = "Library Size", 
          col = SlideName)

plotROIQC(seo_qc_roi, x_threshold = 150, y_threshold = 1e+01, color = SlideName)

# Relative log expression distribution
plotRLExpr(seo) #RLE of raw count 
plotRLExpr(seo_qc_roi)
#Remove the technical variations due to the library size differences
plotRLExpr(seo_qc_roi, ordannots = "SlideName", assay = 2, color = SlideName)
#can also plot by tissue type or other classification
plotRLExpr(seo_qc_roi, ordannots = "tissue", assay = 2, color = tissue)

#Dimentionality Reduction
#seo_qc_roi@assays #Assay 2 is based on logcounts
# BiocManager::install("scater")
drawPCA(seo_qc_roi, assay = 2, color = tissue) # however since the pca will change axis every time we plot, 
#We can save the data to analyze it the same way every time
drawPCA(seo_qc_roi, assay =2, color = SlideName)
#To make it reproducible
set.seed(100)
seoPCA <-  scater::runPCA(seo_qc_roi)
pca_results <-  reducedDim(seoPCA, "PCA")
drawPCA(seoPCA, precomputed = pca_results, col = tissue)
drawPCA(seoPCA, precomputed = pca_results, col = SlideName)

#Draw PCA Scree Plot
plotScreePCA(seo_qc_roi, precomputed = pca_results)
#Plot Pair PCA
plotPairPCA(seo_qc_roi, col= tissue, precomputed = pca_results, n_dimension = 4)
plotPairPCA(seo_qc_roi, col= SlideName, precomputed = pca_results, n_dimension = 4)

#Plot Multidimensional Scaling (MDS) plot
standR::plotMDS(seo_qc_roi, assay = 2, color = SlideName)
standR::plotMDS(seo_qc_roi, assay = 2, color = tissue)


#UMAP
set.seed(100)
seoUMAP <- scater::runUMAP(seoPCA, dimred = "PCA")
plotDR(seoUMAP, dimred = "UMAP", col = tissue)
plotDR(seoUMAP, dimred = "UMAP", col = SlideName)

## Normalization########
#Using the data as seoPCA or seoUMAP produces the same plot, so the data remains same
#As per the vignette, the later is used in the normalization process below.

#names(seoUMAP@metadata)
seo_tmm <- geomxNorm(seoUMAP, method = "TMM") #TMM Normalization
plotRLExpr(seo_tmm, assay = 2, color = tissue) + ggtitle("TMM")
#geomxNorm below, adds two more parameter to the metadata: "norm.factor" "norm.method".
#names(seo_tmm@metadata)
plotRLExpr(seo_tmm, assay = 2, color = SlideName) + ggtitle("TMM")

#Do PCA plot of the geomxNorm'ed data
spe_tmm <- scater::runPCA(seo_tmm)
pca_results_tmm <- reducedDim(spe_tmm, "PCA")
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = SlideName)
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = tissue)

##Batch Correction#########
#Find the least variable 300 genes across the slides
# Dataset used is **from before Normalization** (seoUMAP /seoPCA)
spe_batch <- findNCGs(seoUMAP, batch_name = "SlideName", top_n = 300)
names(spe_batch@metadata) #Adds NCGs (Negative Control Genes to metadata)

#Use RUV4 (remove unwanted variation 4) which requires NCG data to normalize using geomxBatchCorrection
# RUV4 requires : factors of interest, NCGs, Number of unwanted factors to use; smallest number where technical variations no longer exist
# Test out K values until desired level of removal of technical variation
# Optimal K value is that produces a separation of main biology of interest on the PCA plot
# Run Paired PCA plot for k values between 1 to 5

spe_ruv <- geomxBatchCorrection(spe_batch, factors = "tissue", 
                                NCGs = metadata(spe_batch)$NCGs, k =5)
print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = tissue,
                  title = paste0("k = 5 : tissue")))

#Using Tissue as factor in geomBatchCorrection
#Batch Correction 2
#Since Slidename is our biology of interest, batch factor inducing variable could be tissue
#But Tissue is another factor of comparison, and not a batch effect inducer
#So Batch correction is not required

#Trying Alternate Batch Correction Method Limma
#Requires: input object &
#batch: a vector indicating the batch information for all samples;
#design: a design matrix generated from model.matrix, in the design matrix, all biologically-relevant factors should be included.

#Then the different batch correction methods could be compared


##############################
#DEG with limma-voom pipeline (Alternatives: edgeR, DESeq2)
##############################
#Normalized counts should not be used in linear moodelling, rather the weight matrix
#...generated from geomBatchCorrection should be used as covariates
#Weight matrix can be found as following **from Batch Corrected Data**
colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-1, ncol(colData(spe_ruv)))] |>head()
#↑This will have as many W columns as defined in the geomBatchCorrection step

#Establishing design matrix and contrast
#Derive DGElist object from SpatialExperiment Object using SE2DGEList function of edgeR

dge <- SE2DGEList(spe_ruv)
#Check columns and rows of ColData and metadata in the SpatialExperiement Object
colnames(colData(spe_ruv))
names(metadata(spe_ruv))

#Adding W matrices resulting from RUV4 to the model matrix as covariates to use the batch corrected data
design <- model.matrix(~0 + tissue + ruv_W4 + ruv_W5, data = colData(spe_ruv))
colnames(design)
#Rename by removing the tissue prefix
colnames(design) <- gsub("^tissue", "", colnames(design))

#To compare Tumor with Tumor Edge
contr.matrix <- makeContrasts(
  TvE = tumor - tumor_edge,
  levels = colnames(design))

#Reduce the number of genes with low coverage for accurate mean variance relationship
keep <- filterByExpr(dge, design)
table(keep) #Shows no gene is excluded, all 19948 genes kept
rownames(dge)[!keep] #remove any gene with low expression; here nothing will be removed
dge_all <- dge[keep, ]

#Biological CV (BCV) is the coeeficient of variation with with the (unknown) true abundace
#...of the gene varies vetween rerplicate RNA samples.
#BCV Check
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
names(dge_all)

plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))

bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)
 highbcv <- bcv_df$BCV >0.8
 highbcv_df <- bcv_df[highbcv, ]
 points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
 text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)

#Differential Expression
#In the limma-voom pipeline, linear modelling is carried out on the log-CPM values by using 
# ...the voom, lmFit, contrasts.fit and eBayes functions.
v <- voom(dge_all, design, plot = TRUE) #Variance Trend

fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)
summary_efit #Shows that between the TvE contast group, how many genes are up and downregulated

#Visualization
# library(ggrepel)
# library(tidyverse)
de_results_TvE <- topTable(efit, coef = 1, sort.by = "P", n = Inf)

de_genes_toptable_TvE <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

de_results_TvE %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
                     ifelse(logFC<0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
  ggplot(aes(AveExpr, logFC, col = DE)) +
  geom_point(shape = 1, size = 1) +
  geom_text_repel(data = de_genes_toptable_TvE %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
                                       ifelse(logFC < 0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("Tumor vs Tumor_Edge") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5))

#Gene Set Enrichment Analysis
#Using fry from limma package
#Load Gene sets
library(msigdb)
library(GSEABase)
msigdb_hs <- getMsigdb(version ='7.2')

