# Define the project name, organize folder structure. Script in project root.
project <- "nanostring"
# Get the Rstudio Directory
rstudio_dir <- rstudioapi::getActiveProject()
# Define the project directory inside Rstudio directory
project_dir <- file.path(rstudio_dir, project)
if (!file.exists(project_dir)) {dir.create(project_dir, recursive = TRUE)}
# Define input and output directories
input_dir <- file.path(project_dir, "input")
output_dir <- file.path(project_dir, "output")
# Check if input and output directories exist, if not, create them
if (!file.exists(input_dir)) {dir.create(input_dir, recursive = TRUE)}
if (!file.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
# Set the working directory to the project directory
setwd(project_dir)
#Print Confirmation of Correct Folder Structure
if(basename(getwd()) == project) "Folder setup correctly" else "Fix folder structure"
#===
#Vignette using DCC and PKC Files of NanoString: https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
#Vignette using CountData of NanoString: https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
#===
library(tidyverse)
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
#BioProbeCountMmatrix Worksheet: TargetName is default column
featureAnnoFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                                      sheet = "BioProbeCountMatrix") %>% 
  select(., 1:12) %>% 
  select(TargetName, everything())
# Process rows with "NegProbe-WTX"
negprobe_wtx <- featureAnnoFile %>%
  filter(TargetName == "NegProbe-WTX")

# Process other rows to remove duplicates based on "TargetName" and concatenate unique values for duplicates
other_rows <- featureAnnoFile %>%
  filter(TargetName != "NegProbe-WTX") %>%
  mutate_at(c("GeneID"), as.character) %>% 
  group_by(TargetName) %>%
  summarise(across(
    .cols = everything(),
    .fns = ~ if (n() > 1) {
      unique_values <- unique(trimws(unlist(strsplit(as.character(.), ","))))
      paste(unique_values, collapse = ",")
    } else {
      first(.)
    }
  ), .groups = "drop") %>%
  mutate(GeneID = first(featureAnnoFile$GeneID[match(TargetName, featureAnnoFile$TargetName)]))


# Combine both processed datasets
final_featureAnnoFile <- bind_rows(negprobe_wtx, other_rows) %>% 
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
                        sheet="BioProbeCountMatrix") %>% 
  select(., c(3, 13:187)) %>% 
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
assays(seo)
seo$ROILabel
assayNames(seo)

#View the count table using the assay function
assay(seo, "counts")[1:5, 1:5]
assay(seo, "logcounts")[1:5, 1:5]
#Sample metadata is in the colData object
colData(seo)[1:5,1:5]
#Gene metadata stored in the rowData of the object
rowData(seo)[1:5,1:5]
#WTA has NegProbe-WTX data. Ensure that there are no duplicate gene names in the TargetName column
metadata(seo)$NegProbes[,1:5]
colData(seo)$QCFlags

#QC Steps 
#Sample level QC
library(ggplot2)
library(ggalluvial)
plotSampleInfo(seo, column2plot =c("SlideName", "tissue"))#, "CD45", "Neuron"))
# colData(seo)$tissue


#Gene level QC
seo #Dim 19962x175
seo_qc <- addPerROIQC(seo, 
                      sample_fraction = 0.9, #Default
                      rm_genes =TRUE,
                      min_count = 5) #Default
seo_qc #Gene with low count in more than threshold (sample_fraction=0.9)
# are removed by applying the function. Dim 19948x175, so removed 14 genes
dim(seo) #The dimensions remains same even if the row name goes down from 19962 > 19948 
metadata(seo_qc) |> names()
metadata(seo) |>names()
# plotGeneQC(seo, ordannots = "regions", col = regions, point_size = 2)
plotGeneQC(seo_qc, ordannots = "tissue", col = tissue, point_size = 2)
#colData(seo)$regions
#colnames(seo) %>% print() 
data("seo_qc")
seo_subset <-  addPerROIQC(seo_qc)
plotGeneQC(seo_subset)


#ROI level QC
plotROIQC(seo_qc, x_threshold = 150, color = SlideName)
colData(seo_qc)$AOINucleiCount
#AOIINuclei count of 150 looks like a good threshold from the figure
qc <- colData(seo_qc)$AOINucleiCount > 150
table(qc) # 3 Values Below threshold
seo_qc # Dim 19948x175
seo_qc_roi <- seo_qc[, qc]
seo_qc_roi # We removed 3 ROI froom dataset. Dim 19948x172
# Comparing the Library Size with ROI Area size
plotROIQC(seo_qc_roi, 
          x_threshold = 20000, 
          x_axis = "AOISurfaceArea", 
          x_lab = "AreaSize", 
          y_axis = "lib_size", 
          y_lab = "Library Size", 
          col = SlideName)
plotROIQC(seo_qc_roi, x_threshold = 150, color = SlideName)
# Relative log expression distribution
plotRLExpr(seo) #RLE of raw count 
plotRLExpr(seo_qc_roi)
#Remove the technical variations due to the library size differences
plotRLExpr(seo_qc_roi, ordannots = "SlideName", assay = 2, color = SlideName)
#can also plot by tissue type or other classification
plotRLExpr(seo_qc_roi, ordannots = "tissue", assay = 2, color = tissue)

#Dimentionality Reduction
# BiocManager::install("scater")
drawPCA(seo_qc_roi, assay = 2, color = tissue) # however since the pca will change axis everytime we plot we can save the data to analyze it the same way everytime
#To make it reproducible
set.seed(100)
seoPCA <-  scater::runPCA(seo_qc_roi)
pca_results <-  reducedDim(seoPCA, "PCA")
drawPCA(seoPCA, precomputed = pca_results, col = tissue)
drawPCA(seoPCA, precomputed = pca_results, col = SlideName)
#UMAP
set.seed(100)
seoUMAP <- scater::runUMAP(seoPCA, dimred = "PCA")
plotDR(seoUMAP, dimred = "UMAP", col = tissue)
plotDR(seoUMAP, dimred = "UMAP", col = SlideName)

## Normalization
seo_tmm <- geomxNorm(seoUMAP, method = "TMM") #TMM Normalization
plotRLExpr(seo_tmm, assay = 2, color = tissue) + ggtitle("TMM")


