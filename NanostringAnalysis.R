#===
#Vignette using DCC and PKC Files of NanoString: https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
#Vignette using CountData of NanoString: https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
#Vignette partial upto data preparation: https://davislaboratory.github.io/standR/articles/standR_introduction.html
#===

##### Define Functions########################
setupProject <- function(project) {
  #Create 'project' dir if not same as name of *.Rproj dir as root/proj/io, else
  #..set root as proj dir, create io as root/io; set input_dir & output_dir vars
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
install_and_load_packages <- function(cran_packages, bioc_packages) {
  # Install missing CRAN packages
  new_packages_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
  if (length(new_packages_cran) > 0) {install.packages(new_packages_cran)}
  # Install missing Bioconductor packages
  new_packages_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
  if (length(new_packages_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install(new_packages_bioc, update = FALSE)
  }
  # Load all packages
  all_packages <- c(cran_packages, bioc_packages)
  sapply(all_packages, require, character.only = TRUE)
}
sentence_case <- function(name) { # Sentence case first word if not uppercase with/out numbers/"-" (eg.DN-A1)
  # Split the sentence into words
  words <- unlist(strsplit(name, " "))
  # Check if the first word should be converted
  first_word <- words[1]
  if (!grepl("^[A-Z0-9]+$", first_word) && !grepl("-", first_word)) {
    # Convert the first word to sentence case
    words[1] <- paste0(toupper(substring(first_word, 1, 1)), tolower(substring(first_word, 2)))
  }
  # Join the words back into a sentence
  return(paste(words, collapse = " "))
}
heatmap_func_ss <- function(Input_DF,COL_ANNOT, NoofRows){
  Heatmap_Data <- Input_DF
  Heatmap_Data <- Heatmap_Data[order(matrixStats:::rowVars(as.matrix(Heatmap_Data)),decreasing = T),]
  Heatmap_Data <- Heatmap_Data[c(1:No_Of_Rows),]
  Heatmap_Data <- Heatmap_Data[rowSums(is.na(Heatmap_Data)) != ncol(Heatmap_Data), ]
  Heatmap_Data <- as.data.frame(t(scale(t(Heatmap_Data))))
  Heatmap_Data <- as.data.frame(t(Heatmap_Data))
  col = 1
  for(col in length(Heatmap_Data[1,]):1){
    print(col)
    Column_Name = colnames(Heatmap_Data)[col]
    print(Column_Name)
    Heatmap_Data <- Heatmap_Data[order(Heatmap_Data[,Column_Name],decreasing = T),]
  }
  
  Heatmap_Data <- as.data.frame(t(Heatmap_Data))
  
  col = 1
  for(col in length(Heatmap_Data[1,]):1){
    print(col)
    Column_Name = colnames(Heatmap_Data)[col]
    print(Column_Name)
    Heatmap_Data <- Heatmap_Data[order(Heatmap_Data[,Column_Name],decreasing = T),]
  }
  
  
  ##########################
  Heatmap_Annotation_Data <- COL_ANNOT
  ################################
  REQD_SAMPLES <- as.character(intersect(rownames(COL_ANNOT),colnames(Heatmap_Data)))
  ################################
  Heatmap_Annotation_Data <- Heatmap_Annotation_Data[REQD_SAMPLES,]
  Heatmap_Data <- Heatmap_Data[,REQD_SAMPLES]
  ##########################3
  
  Color_Sets <- list(brewer.pal(12,"Set1"),
                     brewer.pal(12,"Accent"),
                     brewer.pal(12,"Paired"),
                     brewer.pal(12,"Dark2")
  )
  
  
  Color_Set1 <- Color_Sets[[1]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"ScanLabel"]),"ScanLabel"])))]
  names(Color_Set1) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"ScanLabel"]),"ScanLabel"])
  
  Color_Set2 <- Color_Sets[[2]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Neuron"]),"Neuron"])))]
  names(Color_Set2) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Neuron"]),"Neuron"])
  
  Color_Set3 <- Color_Sets[[3]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"tissue"]),"tissue"])))]
  names(Color_Set3) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"tissue"]),"tissue"])
  
  Color_Set4 <- Color_Sets[[4]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"CD45"]),"CD45"])))]
  names(Color_Set4) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"CD45"]),"CD45"])
  
  
  
  HEATMAP_ANNOTAION <- HeatmapAnnotation(ScanLabel=as.character(Heatmap_Annotation_Data[,"ScanLabel"]),
                                         Neuron=as.character(Heatmap_Annotation_Data[,"Neuron"]),
                                         tissue=as.character(Heatmap_Annotation_Data[,"tissue"]),
                                         CD45=as.character(Heatmap_Annotation_Data[,"CD45"]),
                                         col = list(ScanLabel=Color_Set1,
                                                    Neuron=Color_Set2,
                                                    tissue=Color_Set3,
                                                    CD45=Color_Set4))
  
  
  
  
  col_fun = colorRamp2(c(min(Heatmap_Data,na.rm = T), 0, max(Heatmap_Data,na.rm = T)), c("blue", "white", "red"))
  
  Heatmap_Data <- Heatmap_Data[,rownames(Heatmap_Annotation_Data)]
  
  INPUT_DATA_HEATMAP <- Heatmap(as.matrix(Heatmap_Data),
                                #width = unit(40, "cm"),height = unit(40, "cm"),
                                #right_annotation = Heatmap_rowAnnotation,
                                top_annotation = HEATMAP_ANNOTAION,    
                                border_gp = gpar(col = "black", lwd = 2),
                                grid.rect(gp = gpar(lwd = 2, fill = "transparent")),
                                name = UNIT,
                                #row_split = IC50_HEATMAP_ROWSPLIT[,"RowOrder"],
                                column_split = Heatmap_Annotation_Data[,c("Group")],
                                rect_gp = gpar(col = "black", lwd =0.2),
                                row_gap = unit(5, "mm"),
                                column_title = " ",
                                column_title_side = "top",
                                row_title = " ", row_title_rot = 90,
                                col = col_fun,
                                row_names_gp = gpar(fontsize = 12,just = "center",fontface ="bold"),
                                column_names_gp = gpar(fontsize = 14,just = "center",fontface ="bold"),
                                column_title_gp = gpar(fontsize = 14,just = "center",fontface ="bold"),
                                row_title_gp = gpar(fontsize = 14,just = "center",fontface ="bold"),
                                column_names_side = "bottom",
                                na_col = "grey90",
                                cluster_rows = T,
                                cluster_columns = T,
                                show_row_dend = T,
                                show_column_dend = T,
                                column_names_centered = TRUE,
                                column_names_max_height = unit(12, "cm"),
                                show_row_names = T,
                                row_names_side = "left",
                                column_names_rot = 90,
                                show_column_names = T,
                                heatmap_legend_param = list(direction = "vertical")
  )
  
  
  return(INPUT_DATA_HEATMAP)
}

##### Setup Project ----
## Initiate project
setupProject("NanostringAnalysis") ; print(paste0("Working dir is: ", getwd()))
# If any project specific override: output folder if any
# output_dir <- paste0(output_dir, "/1.1_Nanostring")

## Install & Load Packages
cran_packages <- c("circlize", "colorRamp2", "DT", "ggalluvial", "ggrepel", "grid", "igraph", "magick", "patchwork", "RColorBrewer", "tidyverse")
bioc_packages <- c("ComplexHeatmap","edgeR", "fsgsea", "GSEABase", "GSVA", "limma", "msigdb", "msigdbr", "qusage", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
install_and_load_packages(cran_packages, bioc_packages)

#### Source & Process Input files ----
# sampleAnnoFile from SegmentProperties sheet
# countFile & featureAnnoFile from BioProbeCountMatrix sheet
#sampleAnnoFile----
#SegmentProperties Worksheet: SegmentDisplayName is default column
sampleAnnoFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                                     sheet="SegmentProperties")

final_sampleAnnoFile <- sampleAnnoFile %>% 
  mutate(SlideName = gsub(" breast", "", #Format texts; remove extra words and signs
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

featureAnnoFile <- dplyr::select(featureAnnoFile, 1:12) %>% #Get the featureData relevant columns only
  dplyr::select(TargetName, everything())

#Process the duplicates where the Targetname is not NegProbe-WTX and then combine them back
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

#==
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

#### countFile----
countFile <- readxl::read_excel("input/Initial Dataset 5-9-23.xlsx", 
                        sheet="BioProbeCountMatrix") 
countFile <- countFile %>%   dplyr::select(., c(3, 13:187)) %>% #Select the rest of sheet after the featureAnnoFile
  # filter(TargetName !=duplicates) %>%
  setNames(c(gsub(" 1/2-", "", colnames(.)))) %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#Process the duplicates where the Targetname is not NegProbe-WTX and then combine them back
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

#Old Code Skip----
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

####### CREATE SPATIAL EXPERIMENT OBJECT----
# BiocManager::install("standR")
# library(standR)
#install.packages("ggalluvial") #required for the standR package
#install.packages("magick") #required for the standR package

seo <- readGeoMx(countFile = final_countFile,
                  sampleAnnoFile = final_sampleAnnoFile,
                  featureAnnoFile = final_featureAnnoFile)

# seo_df <- readGeoMx(countFile = "output/countFile.txt",
                 # sampleAnnoFile = "output/sampleAnnoFile.txt",
                 # featureAnnoFile = "output/featureAnnoFile.txt",
                 # # colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"), #optional as we are using standard
                 # rmNegProbe = TRUE,
                 # NegProbeName = "NegProbe-WTX")

# library(SpatialExperiment)
seo

###Examine the data (Optional: To understand the structure)----
assays(seo)
seo$SlideName
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
names(colData(seo))
#Gene metadata stored in the rowData of the object
rowData(seo)[1:5,1:5]
#WTA has NegProbe-WTX data. Ensure that there are no duplicate gene names in the TargetName column
metadata(seo)$NegProbes[,1:5]
colData(seo)$QCFlags
colData(seo)$tissue
names(seo@colData)

####### QC Steps ---- 
#Sample level QC
# library(ggplot2)
# library(ggalluvial)
#Visualize the data
plotSampleInfo(seo, column2plot =c("SlideName", "tissue", "CD45", "Neuron"))
plotSampleInfo(seo, column2plot =c("SlideName", "tissue")) + coord_flip()
#### Gene level QC----
seo #Dim 19962x175
names(colData(seo)) #46 Columns in ColData
seo_qc <- addPerROIQC(seo, 
                      rm_genes =TRUE,
                      sample_fraction = 0.9, #Default
                      min_count = 5) #Default
seo_qc #Gene with low count and expression values in more than threshold (sample_fraction=0.9)
# are removed by applying the function. Dim 19948x175, so removed 14 genes
dim(seo) 
dim(seo_qc) #Genes not meeting the above criteria were removed 19962 > 19948 
names(seo_qc@colData)
seo_qc@colData
view(colData(seo_qc)[ , c("countOfLowEprGene", "percentOfLowEprGene", "ScanLabel", "lib_size", "tissue")])
view(colData(seo_qc))

length(seo@metadata$genes_rm_rawCount) #0
length(seo_qc@metadata$genes_rm_rawCount) #175
names(seo_qc@metadata)
metadata(seo_qc) |> names() #Same as above
metadata(seo) |>names()

#addPerROIQC added columns to colData :lib_size, countOfLowEprGene, percentOfLowEprGene
# and also added columns to metadata: lcpm_threshold, genes_rm_rawCount, genes_rm_logCPM  


# plotGeneQC(seo_qc, ordannots = "regions", col = regions, point_size = 2)
plotGeneQC(seo_qc)
plotGeneQC(seo_qc, top_n=12, ordannots = "SlideName", col = SlideName, point_size = 2)
plotGeneQC(seo_qc, top_n=12, ordannots = "tissue", col = tissue, point_size = 2)
sapply(seo_qc@colData, class)


#colData(seo)$regions
#colnames(seo) %>% print() 

# data("seo_qc")
# seo_subset <-  addPerROIQC(seo_qc) #Wrong; already done this step to get seo_qc
# plotGeneQC(seo_subset)

#Final data of this segment: seo_qc
#### ROI level QC ----
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

plotROIQC(seo_qc_roi, 
          x_threshold = 150, 
          x_axis = "AOINucleiCount", #Default
          y_threshold = 1e+01, 
          color = SlideName)

# Relative log expression distribution
plotRLExpr(seo) #RLE of raw count 
plotRLExpr(seo_qc_roi)
#Remove the technical variations due to the library size differences
plotRLExpr(seo_qc_roi, ordannots = "SlideName", assay = 2, color = SlideName)
#can also plot by tissue type or other classification
plotRLExpr(seo_qc_roi, ordannots = "tissue", assay = 2, color = tissue)

#### Dimentionality Reduction ----
#seo_qc_roi@assays #Assay 2 is based on logcounts
# BiocManager::install("scater")
drawPCA(seo_qc_roi, assay = 2, color = tissue) # however since the pca will change axis every time we plot, 
#We can save the data to analyze it the same way every time
drawPCA(seo_qc_roi, assay =2, color = SlideName)
#To make it reproducible
set.seed(100)

seoPCA <-  scater::runPCA(seo_qc_roi)
#runPCA adds reducedDimNames PCA
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
#runUMAP adds one more reducedDimNames UMAP
plotDR(seoUMAP, dimred = "UMAP", col = tissue)
plotDR(seoUMAP, dimred = "UMAP", col = SlideName)

names(seoUMAP@metadata)

# ↓ ↓ Q3 Normalization and Saving Data~~~~~~~~~~~~~~~~~~~ ----
# Input: seoUMAP; output: seoUMAP_Normalized_Q3
Counts_Data <- as.data.frame(counts(seoUMAP))
##(SKIP IF NOT SAVING Count Data)Generate a Gene Column for Data Saving ----
Counts_Data <- as.data.frame(cbind(Gene=rownames(Counts_Data), Counts_Data))
write.table(Counts_Data, paste(output_dir, "_Counts.txt", sep=""), 
            sep="\t", row.names = F, col.names = T, quote = T)
#Read Data
Counts_Data <- read.delim(paste(output_dir, "_Counts.txt", sep=""), 
                          sep="\t", header = T, check.names = F, stringsAsFactors = F)
#Assign the rownames back as Genes
rownames(Counts_Data) <- Counts_Data[,"Gene"]
#Remove the Gene Column to make original counts
Counts_Data[,"Gene"] <- NULL

### Add Normalized Data to logcounts field----
NORMALIZATION <- "upperquartile"

seoUMAP_Normalized_Q3 <- geomxNorm(seoUMAP, method=NORMALIZATION, log = T)
names(seoUMAP_Normalized_Q3@metadata) #Added "norm.method" and "norm.factor" to metadata
seoUMAP_Normalized_Q3@metadata$norm.method

plotRLExpr(seoUMAP_Normalized_Q3, assay = 2, color = tissue) + ggtitle("Q3 Normalized")
plotRLExpr(seoUMAP_Normalized_Q3, assay = 2, color = SlideName) + ggtitle("Q3 Normalized")

Normalized_Counts_Data_Q3 <- as.data.frame(logcounts(seoUMAP_Normalized_Q3))

##(SKIP IF NOT SAVING Normalized Count Data)Generate a Gene Column for Data Saving ----
Normalized_Counts_Data <- as.data.frame(cbind(Gene=rownames(Normalized_Counts_Data), Normalized_Counts_Data))
write.table(Normalized_Counts_Data, paste(output_dir, NORMALIZATION, "_Normalized_Counts.txt", sep=""), 
            sep="\t", row.names = F, col.names = T, quote = T)
#Read Data
Normalized_Counts_Data <- read.delim(paste(output_dir, NORMALIZATION, "_Normalized_Counts.txt", sep=""), 
            sep="\t", header = T, check.names = F, stringsAsFactors = F)
#Assign the rownames back as Genes
rownames(Normalized_Counts_Data) <- Normalized_Counts_Data[,"Gene"]

#Remove the Gene Column to make original normalized counts
Normalized_Counts_Data[,"Gene"] <- NULL

##Make Feature Data ----
Feature_Data <- as.data.frame(rowData(seoUMAP_Normalized_Q3))
#Skip if not saving
# Feature_Data <- as.data.frame(cbind(Gene=rownames(Feature_Data), Feature_Data))

Sample_Data <- as.data.frame(colData(seoUMAP_Normalized_Q3))
#Skip if not saving
#Sample_Data <- as.data.frame(cbind(ROI=rownames(Sample_Data), Sample_Data))


#↑ ↑ Q3 Normalization End~~~~~~~~~~~~~~~~~~

#### TMM Normalization########----
# Input: seoUMAP output: seoUMAP_Normalized_TMM
#Data has become seoUMAP now with addition of two reducedDimNames PCA and UMAP
#names(seoUMAP@metadata)
NORMALIZATION_METHOD <- "TMM"
seoUMAP_Normalized_TMM <- geomxNorm(seoUMAP, method = NORMALIZATION_METHOD) #TMM Normalization
#geomxNorm Adds two more items to metadata: norm.factor, norm.method

#names(seoUMAP_Normalized_TMM@metadata)
plotRLExpr(seoUMAP_Normalized_TMM, assay = 2, color = tissue) + ggtitle("TMM")
plotRLExpr(seoUMAP_Normalized_TMM, assay = 2, color = SlideName) + ggtitle("TMM")
Normalized_Counts_Data_TMM <- as.data.frame(logcounts(seoUMAP_Normalized_TMM))
##↓↓~~~~~~~~~~~~~~~~~~~~~
####Plot QC Figures ----
# library(patchwork)

ROI_ANNOTATION_COLS <- c("SlideName", "ScanLabel", "Neuron", "tissue", "CD45")
dictionary <- c(TMM=seoUMAP_Normalized_TMM, Q3=seoUMAP_Normalized_Q3) #seoUMAP_Normalized_TMM (for TMM) or seoUMAP_Normalized_Q3 (for Q3)
# eval(GRP)
i = 1
NORM = "Q3" #TMM or Q3 for naming the plots and files

for(i in 1:length(ROI_ANNOTATION_COLS)){
  GRP <- ROI_ANNOTATION_COLS[[i]]
  print(GRP)
  set.seed(100)
  RLE_PLOT <- plotRLExpr(dictionary[[NORM]], assay = 2, ordannots=GRP, color=get(GRP)) + ggtitle(paste(NORM, GRP, "RLE")) + labs(color=GRP)
  PCA_SCREE_PLOT <- plotScreePCA(dictionary[[NORM]], assay=2, dims=10)
  PAIRWISE_PCA <- plotPairPCA(dictionary[[NORM]], title=GRP, col=get(GRP), shape=get(GRP), assay=2, n_dimension=2) + labs(color=GRP, shape=GRP)
  PCA_BI_PLOT <- plotPCAbiplot(dictionary[[NORM]], n_loadings=10, assay=2, col=get(GRP)) + labs(color=GRP)
  MultiDIM_SCALING_PLOT <- standR::plotMDS(dictionary[[NORM]], assay=2, col=get(GRP), shape=get(GRP)) + labs(color=GRP, shape=GRP)
  UMAP <- plotDR(dictionary[[NORM]], dimred="UMAP", col=get(GRP)) + labs(color=GRP)
  
  graphics.off() #Set up the PDF output
  pdf(paste(output_dir,"/", NORM, "_NORMALIZED_", GRP, ".pdf", sep = ""), width=15, height=20)
  combined_plot <- (
    (RLE_PLOT / PAIRWISE_PCA) /
      (PCA_BI_PLOT + MultiDIM_SCALING_PLOT + UMAP) / PCA_SCREE_PLOT
  ) +
    plot_layout(ncol = 1, guides = 'collect')
  
  print(combined_plot) #Print the combined plot to the PDF
  graphics.off() # Close the PDF device
}

#### ssGSEA Analysis----
# library("RColorBrewer")
# install.packages("colorRamp2")
# library("colorRamp2")
# BiocManager::install("qusage")
# library("qusage")

MSIG_DB <- paste0(input_dir, "/MSIG_DB/")
GeneSets <- list.files(MSIG_DB, pattern= "*.gmt", full.names = F)
MSigDB_Dictionary <- list(mh_all="HallMark Gene Sets", 
                          m2_cp_biocarta="BioCarta subset of Canonical Pathways",
                          m2_cp_reactome="Reactome subset of Canonical Pathways",
                          m2_cp_wikipathways="WikiPathways subset of Canonical Pathways",
                          m8_all="Cell Type Signature Gene Sets")
#MSigDB_Dictionary$m8.all

geneset =1
for (geneset in 1:length(GeneSets)){
  # geneset_name <-  "m2_cp_biocarta" #For Non-loop version
  geneset_name = gsub("\\.","_",
                      gsub(".v2023.2.Mm.symbols.gmt", "", GeneSets[geneset]))
  
  print(geneset_name)
  ####
  gset_mouse=paste(MSIG_DB, GeneSets[geneset], sep = "")
  # Signature <- qusage::read.gmt(paste(MSIG_DB,"m2.cp.biocarta.v2023.2.Mm.symbols.gmt", sep="")) #Non-loop version
  Signature <- qusage::read.gmt(gset_mouse)
  ####
  
  ssParam <- gsvaParam(as.matrix(Normalized_Counts_Data), #Matrix of gene expression
                       Signature, #List object of gene sets
                       kcdf = "Gaussian", #Default
                       maxDiff = TRUE) #Default
  
  ssGSEAScores <- GSVA::gsva(ssParam, verbose=TRUE)
  dim(ssGSEAScores)
  
  # ssGSEAScores <- gsva(as.matrix(Normalized_Counts_Data), #Edited for the new version of gsva package
  #                      Signature, 
  #                      method="ssgsea",
  #                      ssgsea.norm=FALSE,
  #                      kcdf="Gaussian")
  ssGSEAScores <- as.data.frame(ssGSEAScores)
  
  ##(SKIP IF NOT SAVING ssGSEAScores.txt)
  # ssGSEAScores <- as.data.frame(cbind(GSET=rownames(ssGSEAScores), ssGSEAScores))
  # write.table(ssGSEAScores, paste(output_dir, geneset_name, "_ssGSEAScores.txt", sep=""),
  #             sep = "\t", row.names = F, col.names = T)
  # INPUT_DF <- NULL
  # INPUT_DF <-  ssGSEAScores
  # rownames(INPUT_DF) <- INPUT_DF[,1]
  # INPUT_DF[,1] <- NULL

  Input_DF <- ssGSEAScores
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("HALLMARK_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("REACTOME_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("BIOCARTA_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("WP_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("TABULA_MURIS_SENIS_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("DESCARTES_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("ZHANG_","",as.character(x)))
  ###########################
  COL_ANNOT <- NULL
  COL_ANNOT = Sample_Data[,c("ROI","ScanLabel","Neuron","tissue","CD45")]
  COL_ANNOT[,"Group"] <- COL_ANNOT[,"ScanLabel"]
  rownames(COL_ANNOT) <- COL_ANNOT[,"ROI"]
  COL_ANNOT[,"ROI"] <- NULL
  ##########################
  No_Of_Rows = 100
  # UNIT = geneset_name
  UNIT = get("MSigDB_Dictionary")[geneset_name]
  #########################3

  ssGSEA_Heatmap <- heatmap_func_ss(Input_DF,COL_ANNOT,100)
  
  #NCOLS = as.character(length(Heatmap_Data[,1]))
  #NROWS = as.character(length(Heatmap_Data[1,]))
  #######################################
  graphics.off()
  pdf(paste(output_dir, "/", geneset_name,"_ssGSEAScores.pdf",sep=""),width =60,height =30)
  draw(ssGSEA_Heatmap,padding = unit(c(1,5,1,1), "in"),heatmap_legend_side = "right",row_title = "", row_title_gp = gpar(col = "red"),legend_grouping = "original",
       column_title = paste("ssGSEA- ", UNIT,sep=""), column_title_gp = gpar(fontsize = 32))
  graphics.off()
  #####################################
}



##↑↑~~~~~~~~~~~~~~~~~~~~~~ ----

### ↓ MY Pathway specific codes start. *****************************



### ↑ MY Pathway specific codes end. *****************************

#Do PCA plot of the geomxNorm'ed data (to see only, not intended to add data to object)
seoUMAP_Normalized_TMM <- scater::runPCA(seoUMAP_Normalized_TMM) #Adds pca data to reducedDim
pca_results_tmm <- reducedDim(seoUMAP_Normalized_TMM, "PCA")
plotPairPCA(seoUMAP_Normalized_TMM, precomputed = pca_results_tmm, color = SlideName)
plotPairPCA(seoUMAP_Normalized_TMM, precomputed = pca_results_tmm, color = tissue)

#### Batch Correction########----
#### Input: seoUMAP_Normalized_TMM/seoUMAP_Normalized_Q3 -> Output: seoNCG_TMM/seoNCG_Q3
# Dataset used is **from before Normalization** (seoUMAP that has two reduced DimNames added)
# One way to remove batch effect is to identify negative control genes (NCGs) among the slides
# Find the least variable 300 genes across the slides; using findNCGs function
seoNCG_TMM <- findNCGs(seoUMAP_Normalized_TMM, 
                      n_assay = 2, #Of two assays in seoUMAP@assays, defining 1st(counts) or 2nd (logcounts) to use
                      batch_name = "SlideName", #Where batch effect expected to be observed
                      top_n = 300)
##findNCGs adds 'NCGs' (Negative Control Genes to metadata)
names(seoNCG_TMM@metadata)
names(seoNCG_TMM@assays)
seoUMAP@colData$Neuron
seoUMAP@colData$SlideName
seoUMAP@colData$tissue

# Use RUV4 (remove unwanted variation 4) which requires NCG data to normalize using geomxBatchCorrection
# RUV4 requires : factors of interest, NCGs, Number of unwanted factors to use; smallest number where technical variations no longer exist
# Test out K values until desired level of removal of technical variation
# Optimal K value is that produces a separation of main biology of interest on the PCA plot
# Run Paired PCA plot for k values between 1 to 5

####Determine parameter to batch correct----
# Input seoNCG_TMM -> seoBatch_TMM -> spe_TMM
#To plot all pairPCAs together
for(i in seq(5)){
  seoBatch_TMM <- geomxBatchCorrection(seoNCG_TMM, factors = "Neuron", #Or use other variables like tissue, Neuron, SlideName
                                  NCGs = metadata(seoNCG_TMM)$NCGs, k = i)
  
  print(plotPairPCA(seoBatch_TMM, assay = 2, n_dimension = 4, color = Neuron, title = paste0("k = ", i)))
  
}

#Else, print each one separately by changing k each time
seoBatch_TMM <- geomxBatchCorrection(seoNCG_TMM, factors = "Neuron", #factors of biological interest
                                NCGs = metadata(seoNCG_TMM)$NCGs, k =2)
## geomxBatchCorrection adds colData 5 data points: ruv_W1 to ruv_W5

print(plotPairPCA(seoBatch_TMM, assay = 2, n_dimension = 4, color = Neuron,
                  title = paste0("k = 2 : Neuron"))) #Change k in plot for each loop manually

#Identify which k has the best separation of the expected biological variation. Here we decide k=4

####Prepare batch corrected data
#Using that k value make the seoBatch_TMM object final, and replace the PCA data with a new PCA
seoBatch_TMM <- geomxBatchCorrection(seoNCG_TMM, factors = "Neuron",  #Select any factor of biological interest such as tissue or Neuron
                                NCGs = metadata(seoNCG_TMM)$NCGs, k = 2)
set.seed(100)
seoBatch_TMM <- scater::runPCA(seoBatch_TMM) #update seoBatch_TMM with new PCA data

pca_results_ruv <- reducedDim(seoBatch_TMM, "PCA") #Extract PCA data
plotPairPCA(seoBatch_TMM, precomputed = pca_results_ruv, color=tissue,
            title = "RUV, k=4")
spe_TMM <- seoBatch_TMM
#Using Tissue as factor in geomBatchCorrection instead of SlideNames in the findNCGs()
#May be alternate approach, however
#Since SlideNames is our biology of interest, batch factor inducing variable could be tissue
#But Tissue is another factor of comparison, and not a batch effect inducer
#So I feel, Batch correction is not required

#Trying Alternate Batch Correction Method Limma
#Requires: input object &
#batch: a vector indicating the batch information for all samples;
#design: a design matrix generated from model.matrix, 
#in the design matrix, all biologically-relevant factors should be included.

#Then the different batch correction methods could be compared


##############################
#DEG with limma-voom pipeline (Alternatives: edgeR, DESeq2)
##############################
#Normalized counts should not be used in linear modelling, rather the weight matrix
#...generated from geomBatchCorrection should be used as covariates
#Weight matrix can be found as following **from Batch Corrected Data**
colData(spe_TMM)[,seq(ncol(colData(spe_TMM))-1, ncol(colData(spe_TMM)))] |>head()
#↑This will have as many W columns as defined in the geomBatchCorrection step

#spe_TMM@colData
temp <-  spe_TMM@colData
spe_TMM@colData$temp

#Establishing design matrix and contrast
#Derive DGElist object from SpatialExperiment Object using SE2DGEList function of edgeR

dge <- SE2DGEList(spe_TMM)
#Check columns and rows of ColData and metadata in the SpatialExperiement Object
colnames(colData(spe_TMM))
names(metadata(spe_TMM))


##### 1 EdgeR : Workflow Beginning for Each New Set of Design Compare----
#Adding W matrices resulting from RUV4 to the model matrix as covariates to use the batch corrected data
design <- model.matrix(~0 + Neuron + ruv_W1 + ruv_W2, data = colData(spe_TMM)) #Select Any factor for comparison like tissue, Neuron, or SlideName
colnames(design)
#Rename by removing the tissue prefix
# colnames(design) <- gsub("^tissue", "", colnames(design))
colnames(design) <- gsub("^Neuron", "", colnames(design))
colnames(design) <- gsub("-H\\+", "H_POS", colnames(design))
colnames(design) <- gsub("-H\\-", "H_NEG", colnames(design))

#Create alternative Design file for DESeq2
design_mods <- as.data.frame(design) %>%
  rownames_to_column("RowName") %>%        #Convert row names to a column
  mutate(Sample = gsub(" \\|.*", "", RowName)) %>%  #Extract the part before the first pipe
  pivot_longer(cols = tumor:tumor_gland, names_to = "Tissue", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(-Value) %>%
  arrange(RowName) %>%
  column_to_rownames("RowName") #Convert back to row names

#The comparison matrix is setup using makeContrasts function from Limma
#To compare Tumor with Tumor Edge
contr.matrix <- makeContrasts(
  HvN = NFH_NEG - NFH_POS,
  levels = colnames(design))

#Filter out the number of genes with low coverage for accurate mean variance relationship
#..using filterByExpr function of egdeR package
keep <- filterByExpr(dge, design)
table(keep) #Shows 1128 genes are excluded, from the total 19948 genes


rownames(dge)[!keep] #remove any gene with low expression
dge_all <- dge[keep, ]
names(dge_all) #Had 3 cols: "counts"  "samples" "genes" 
###QC: Biological CV (BCV) is the coefficient of variation with which the (unknown) true abundance
#...of the gene varies between replicate RNA samples.
#BCV Check
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
names(dge_all) # added 9 cols like "design", "common.dispersion", "trended.dispersion", "tagwise.dispersion", "AveLogCPM", "trend.method", "prior.df", "prior.n", "span"

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

#### 2 Differential Expression LIMMA VOOM ----
#In the limma-voom pipeline, linear modelling is carried out on the log-CPM values by using 
# ...the voom, lmFit, contrasts.fit and eBayes functions.
v <- voom(dge_all, design, plot = TRUE) #Variance Trend

fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)
summary_efit #Shows that between the HvN contast group, how many genes are up and downregulated

#### Volcano Visualization
# library(ggrepel)
# library(tidyverse)
de_results_HvN <- topTable(efit, coef = 1, sort.by = "P", n = Inf)

de_genes_toptable_HvN <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

de_results_HvN %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
                     ifelse(logFC<0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
  ggplot(aes(AveExpr, logFC, col = DE)) +
  geom_point(shape = 1, size = 1) +
  geom_text_repel(data = de_genes_toptable_HvN %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
                                       ifelse(logFC < 0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("NF-H-Negative vs NF-H-Positive") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~----
#### DUPLICATE Second set of comparison Between SlideNames COMMENTED> MAKE FUNCTION ----



# ##### 2 DUPLICATE SET from Above: BEGIN: SlideName: Workflow Beginning for Each New Set of Design Compare--
# #Adding W matrices resulting from RUV4 to the model matrix as covariates to use the batch corrected data
# design <- model.matrix(~0 + SlideName + ruv_W1 + ruv_W2, data = colData(spe_TMM)) #Select Any factor for comparison like tissue, Neuron
# colnames(design)
# #Rename by removing the tissue prefix
# colnames(design) <- gsub("^SlideName", "", colnames(design))
# # colnames(design) <- gsub("-H\\+", "H_POS", colnames(design))
# # colnames(design) <- gsub("-H\\-", "H_NEG", colnames(design))
# 
# #Create alternative Design file for DESeq2
# design_mods <- as.data.frame(design) %>%
#   rownames_to_column("RowName") %>%        #Convert row names to a column
#   mutate(Sample = gsub(" \\|.*", "", RowName)) %>%  #Extract the part before the first pipe
#   pivot_longer(cols = tumor:tumor_gland, names_to = "Tissue", values_to = "Value") %>%
#   filter(Value == 1) %>%
#   select(-Value) %>%
#   arrange(RowName) %>%
#   column_to_rownames("RowName") #Convert back to row names
# 
# ####For Within one group comparisons start here.
# #The comparison matrix is setup using makeContrasts function from Limma
# #To compare Tumor with Tumor Edge
# factor1 = "BBP"
# factor2= "FP"
# 
# title <- paste(factor1, "vs", factor2)  
# # # resultsNames(dds)
# # e <- as.formula(paste0(factor1, "-", factor2))
# 
# 
# contr.matrix2 <- makeContrasts(
#   BBPvsFP = BBP - FP,
#   levels = colnames(design))
# 
# #Filter out the number of genes with low coverage for accurate mean variance relationship
# #..using filterByExpr function of egdeR package
# keep <- filterByExpr(dge, design)
# table(keep) #Shows 893 genes are excluded, from the total 19948 genes (Same for the same column)
# 
# rownames(dge)[!keep] #remove any gene with low expression
# dge_all <- dge[keep, ]
# names(dge_all) #Had 3 cols: "counts"  "samples" "genes" 
# ###QC: Biological CV (BCV) is the coefficient of variation with which the (unknown) true abundance
# #...of the gene varies between replicate RNA samples.
# #BCV Check
# dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
# names(dge_all) # added 9 cols like "design", "common.dispersion", "trended.dispersion", "tagwise.dispersion", "AveLogCPM", "trend.method", "prior.df", "prior.n", "span"
# 
# plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3), main = title)
# 
# bcv_df <- data.frame(
#   'BCV' = sqrt(dge_all$tagwise.dispersion),
#   'AveLogCPM' = dge_all$AveLogCPM,
#   'gene_id' = rownames(dge_all)
# )
# highbcv <- bcv_df$BCV >0.8
# highbcv_df <- bcv_df[highbcv, ]
# points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
# text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)
# 
#### Differential Expression LIMMA VOOM----
# #In the limma-voom pipeline, linear modelling is carried out on the log-CPM values by using 
# # ...the voom, lmFit, contrasts.fit and eBayes functions.
# v <- voom(dge_all, design, plot = TRUE) #Variance Trend
# 
# fit <- lmFit(v)
# fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix2) #Change Contr.Matrix for each Comparison
# efit <- eBayes(fit_contrast, robust = TRUE)
# 
# results_efit <- decideTests(efit, p.value = 0.05)
# summary_efit <- summary(results_efit)
# summary_efit #Shows that between the HvN contast group, how many genes are up and downregulated

# #### Volcano Visualization
# # library(ggrepel)
# # library(tidyverse)
# de_results_HvN <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
# 
# de_genes_toptable_HvN <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)
# 
# de_results_HvN %>% 
#   mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
#                      ifelse(logFC<0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
#   ggplot(aes(AveExpr, logFC, col = DE)) +
#   geom_point(shape = 1, size = 1) +
#   geom_text_repel(data = de_genes_toptable_HvN %>% 
#                     mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP",
#                                        ifelse(logFC < 0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>% 
#                     rownames_to_column(), aes(label = rowname)) +
#   theme_bw() +
#   xlab("Average log-expression") +
#   ylab("Log-fold-change") +
#   ggtitle(title) + #Change for Each Comparison
#   scale_color_manual(values = c("blue", "grey", "red")) +
#   theme(text = element_text(size=15),
#         plot.title = element_text(hjust = 0.5))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~----
############Gene Set Enrichment Analysis LIMMA----
#Using fry from limma package
#Load Gene sets
# library(msigdb)
# library(GSEABase)
msigdb_hs <- getMsigdb(version ='7.2')

msigdb_hs <- appendKEGG(msigdb_hs)

sc <- listSubCollections(msigdb_hs)
gsc <- c(subsetCollection(msigdb_hs, c('h')),
         subsetCollection(msigdb_hs, 'c2', sc[grepl("^CP:", sc)]),
         subsetCollection(msigdb_hs, 'c5', sc[grepl("^GO:", sc)])) %>%
  GeneSetCollection()

#Enrichment Analysis
fry_indices <- ids2indices(lapply(gsc, geneIds), rownames(v), remove.empty = FALSE)
names(fry_indices) <- sapply(gsc, setName)

gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
gsc_category <- gsc_category[sapply(fry_indices, length) > 5] #Filter out genesets with <5 genes

gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]

fry_indices <- fry_indices[sapply(fry_indices, length) > 5]
names(gsc_category) = names(gsc_subcategory) = names(fry_indices)

#Run fry with all the gene sets filtered above
fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
fry_res_out <- lapply(fry_indices_cat, function (x) {
  limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
})

post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
  names(fry_output) <- NULL
  fry_output <- do.call(rbind, fry_output)
  fry_output$GenesetName <- rownames(fry_output)
  fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
  fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
  return(fry_output)
}

fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
  as.data.frame() %>%
  filter(FDR < 0.05) 

#The output data.frame can be inspected for top N genes
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")

#Plotting the top n downregulated genes
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")

############################## ---
#=====DEG with DESeq2 Pipeline ----
############################## ---
## Install & Load Packages (DESeq Specific Packages)
cran_packages <- c("annotate", "circlize", "devtools", "EnhancedVolcano", "ggpubr", "ggrepel", "matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "viridis")
bioc_packages <- c("apeglm", "clusterProfiler", "ComplexHeatmap", "DESeq2", "DOSE", "enrichplot", "genefilter", "GSVA", "org.Hs.eg.db", "org.Mm.eg.db", "pathview", "xCell")
install_and_load_packages(cran_packages, bioc_packages)

# Input SpatialExperiment Object: spe_TMM
fc <- assays(spe_TMM, 1, withDimnames = FALSE)$counts
fc <- assays(spe_TMM, 2, withDimnames =TRUE, stringsAsFactors=FALSE)$counts
sapply(assay(spe_TMM, 2, withDimnames =TRUE)$counts, class)
sapply(countsData, class)

countData <- as.data.frame(counts(spe_TMM))
str(countData)

colData <- as.data.frame(spe_TMM@colData@listData)
colData <- colData %>% dplyr::select(c("SlideName", "Neuron", "tissue", "CD45"))

colData <- colData %>% mutate(tissue = str_to_sentence(tissue))
colData <- colData %>% rename(Tissue = tissue)
colData[,"Neuron"] <- sapply(colData[,"Neuron"], function(x) gsub("-H\\+", "H_POS", as.character(x)))
colData[,"Neuron"] <- sapply(colData[,"Neuron"], function(x) gsub("-H\\-", "H_NEG", as.character(x)))
colData[,"CD45"] <- sapply(colData[,"CD45"], function(x) gsub("\\-", "_NEG", as.character(x)))
colData[,"CD45"] <- sapply(colData[,"CD45"], function(x) gsub("\\+", "_POS", as.character(x)))

str(colData)
#Convert all character columns into factors.
colData <-  colData %>% mutate_if(is.character, as.factor)
str(colData)

# library("DESeq2")
# Create DESeqDataSet object
ddsObject <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ CD45)
# ddsObject <- DESeqDataSet(spe_TMM, design = ~ CD45)

#Keeping rows that have at least 10 reads for a minimum number of samples
#Minimal number of samples is the smallest group size, eg here 12 of each cellLine
#..or minimal number of samples for which non-zero counts would be considered interesting; 3 replicates
# if (perform_subset_analysis) {
#   smallestGroupSize <- 3
# } else {
#   smallestGroupSize <- 12
# }
## Make the previously made spe_TMM as the dds object
ddsObject <- spe_TMM
smallestGroupSize <- 25 #Approx in one group BBP, Sympa etc
#counts(ddsObject)
#keep <-  rowSums(counts(dds2))>= 10
keep <-  rowSums(counts(ddsObject)>= 10) >= smallestGroupSize

ddsObject_filtered <- ddsObject[keep,]

names(ddsObject_filtered@colData)
typeof(ddsObject_filtered@colData$SlideName)
ddsObject_filtered@colData$SlideName
ddsObject_filtered@colData$SlideName <- as.factor(ddsObject_filtered@colData$SlideName)
#Since we are primarily comparing between different groups such as Neuron Slide etc, 
#so our primary level of comparison is that comparison, need to define reference level
#(this only reorders, since the default comparison is with first in the list)
ddsObject_filtered$SlideName <- relevel(ddsObject_filtered$SlideName, ref="BBP")

ddsObject_filtered$SlideName <- droplevels(ddsObject_filtered$SlideName) #remove the levels (of CD45) 
# ...which do not have samples in the current data set. Here nothing removed

##### Run DESeq2 Analysis ----
dds <- DESeq(ddsObject_filtered)

sapply(spe_TMM@assays@data$counts, class)
spe_TMM@assays@data$counts <- as.numeric(unlist(spe_TMM@assays@data$counts))
dds <- DESeqDataSet(spe_TMM)
##### DDS Apply Transformation ----
# Apply transformation & estimate dispersion trend
# vsd <- vst(dds, blind = FALSE) # VST: Variance Stabilizing Transformation
rld <- rlog(dds, blind=FALSE) # RLT: Regularized Log Transformation (Selected for this analysis)

head(assay(rld), 2)
### Heatmaps ----
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#E32600"))(99)

# Extract Transformed values
#rld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds)
# Using rld transformed data for the analysis [rld <- rlog(dds)]
length(rld)

## creating distance matrix (Not Needed)
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),
                       annotation_col = colData, 
                       col=my_colors,
                       main = paste("Distance Matrix"),#, subsetToAnalyze, "Cell Line"),
                       annotation_legend=TRUE)

#Plot Heatmap of Hierarchical Clustering
rld_mat <- assay(rld) #Extract the transformed matrix from the object
rld_cor <- cor(rld_mat) #Compute pairwise correlation values
head(rld_cor)
heat.colors <- RColorBrewer::brewer.pal(6, "BrBG")
pheatmap(rld_cor, 
         # annotation=colData, 
         color = heat.colors, border_color = NA,
         main = paste("Hierarchical Correlation", subsetToAnalyze, "Cell Line"),
         fontsize = 10, fontsize_row = 10)

## PCA Plot using RLT Transformation #Almost a similar plot
plotPCArld <- plotPCA(rld, intgroup="CD45", returnData=FALSE, ntop=length(rld)) 
plotPCArldData <- plotPCA(rld, intgroup="CD45", returnData=TRUE, ntop=length(rld)) 
pcaRLDAll <- plotPCArld+geom_label_repel(data=plotPCArldData, 
                                         aes(label=name), 
                                         min.segment.length = 0.5)+
  ggtitle(label="PCA Plot of All Genes")+
  labs(caption = "Regularized Log Transformed (RLT) Count Data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.caption = element_text(hjust=1, size = '8', color = 'grey', face = 'italic'))

print(pcaRLDAll)

## PCA Plot using RLD for default top 500 genes
plotPcaRLD500<-plotPCA(rld, intgroup="CD45", returnData=FALSE)
plotPcaRLD500Data <- plotPCA(rld, intgroup="CD45", returnData=TRUE)
pcaRLD500 <- plotPcaRLD500 + geom_label_repel(data=plotPcaRLD500Data, 
                                              aes(label=name),
                                              min.segment.length = 0.5,
                                              max.overlaps = 1)+
  ggtitle(label="PCA Plot of Top 500 Variable Genes")+
  labs(caption = "Regularized Log Transformed (RLT) Count Data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5,  size = 14, face = "bold"),
        plot.caption = element_text(hjust=1, size = '8', color = 'grey', face = 'italic'),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "bottom",
        #strip.background = element_blank(),
        axis.text = element_text(size = 10)
        #panel.border = element_rect(color="black"),
  )
print(pcaRLD500)

### DEG Prepare Data for Plot (INPUT NEEDED- Define Comparison Groups)----
ComparisonColumn <- "CD45"
factor1 <- "CD45_POS" #Choose from 128.10, 128.13 & 130
factor2 <- "CD45_NEG"
title <- paste(factor1, "vs", factor2)  
# resultsNames(dds)
e <- as.character(c(ComparisonColumn, factor1, factor2))
#Shrink dds based on comparison data
res <- lfcShrink(dds, contrast = e, type = "normal") #Shrinked Result (L2FC, padj etc); "normal" algorithm

annotation_colors <- list(
  drug = c("128.10"="#9FD900", 
           "128.13"="#FAA800", 
           "130"="#ff5d8f", 
           "Control"="#6c757d"),
  
  cellLine =c("358"="#006E18", 
              "318"="#832161")
)
icolors <- colorRampPalette(c("blue",
                              "white",
                              "red"))(99)
#Combine shrunk results (with lfc, padj etc with normalized count data)
resdata <- merge(as.data.frame(res), #Get the counts data alongwith the comparison results
                 as.data.frame(assay(rld)), #To get transformed count values
                 #as.data.frame(counts(dds, normalized=FALSE)), #To get original count data
                 by = "row.names",
                 sort = FALSE)

names(resdata)[1] <- "EnsembleID" #Rename the first column (Row.names) as EnsembleID
#Get Bioconductor Annotation Database
sp <- org.Mm.eg.db

#Add a column of Gene translated from EnsembleID
resdata$Gene<- mapIds(sp, keys=resdata$EnsembleID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata$EntrezID<- mapIds(sp, keys=resdata$EnsembleID, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")

resdata <- resdata[order(resdata$padj),] #Order by ascending p-adj significance;low p =significant at top
resdata <- resdata[complete.cases(resdata$padj),] #Keep rows which have p-adj data
resdata <- resdata[complete.cases(resdata$Gene),] #Remove rows that don't have assigned genes
resdata <- resdata[!duplicated(resdata$Gene),] #Remove rows with same gene, keeping significant ones (low padj)
resdata <- resdata[order(resdata$log2FoldChange),] #Now order per log2FoldChange :Ascending
