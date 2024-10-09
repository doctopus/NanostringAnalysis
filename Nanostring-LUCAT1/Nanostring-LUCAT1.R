#====
#Vignette using DCC and PKC Files of NanoString: https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
#Vignette using CountData of NanoString: https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
#Vignette partial upto data preparation: https://davislaboratory.github.io/standR/articles/standR_introduction.html
#====

######## Define Functions########################
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
sentence_case <- function(name) {
  # Sentence case first word if not uppercase with/out numbers/"-" (eg.DN-A1)
  # Split the sentence into words
  words <- unlist(strsplit(name, " "))
  # Check if the first word should be converted
  first_word <- words[1]
  if (!grepl("^[A-Z0-9]+$", first_word) && !grepl("-", first_word)) {
    # Convert the first word to sentence case
    words[1] <- paste0(toupper(substring(first_word, 1, 1)), 
                       tolower(substring(first_word, 2)))
  }
  # Join the words back into a sentence
  return(paste(words, collapse = " "))
}
heatmap_func_ss_old <- function(Input_DF,COL_ANNOT, NoofRows){
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
                     brewer.pal(12,"Accent")#,
                     # brewer.pal(12,"Paired"),
                     # brewer.pal(12,"Dark2")
  )
  
  Color_Set1 <- Color_Sets[[1]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Samples"]),"Samples"])))]
  names(Color_Set1) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Samples"]),"Samples"])
  
  Color_Set2 <- Color_Sets[[2]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Tumor"]),"Tumor"])))]
  names(Color_Set2) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Tumor"]),"Tumor"])
  
  # Color_Set1 <- Color_Sets[[1]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"ScanLabel"]),"ScanLabel"])))]
  # names(Color_Set1) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"ScanLabel"]),"ScanLabel"])
  # 
  # Color_Set2 <- Color_Sets[[2]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Neuron"]),"Neuron"])))]
  # names(Color_Set2) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Neuron"]),"Neuron"])
  # 
  # Color_Set3 <- Color_Sets[[3]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"tissue"]),"tissue"])))]
  # names(Color_Set3) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"tissue"]),"tissue"])
  # 
  # Color_Set4 <- Color_Sets[[4]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"CD45"]),"CD45"])))]
  # names(Color_Set4) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"CD45"]),"CD45"])
  
  
  
  HEATMAP_ANNOTAION <- HeatmapAnnotation(Sample=as.character(Heatmap_Annotation_Data[,"Samples"]),
                                         Tissue=as.character(Heatmap_Annotation_Data[,"Tumor"]),
                                         # tissue=as.character(Heatmap_Annotation_Data[,"tissue"]),
                                         # CD45=as.character(Heatmap_Annotation_Data[,"CD45"]),
                                         col = list(Sample=Color_Set1,
                                                    Tumor=Color_Set2))
                                                    # tissue=Color_Set3,
                                                    # CD45=Color_Set4))
  
  
  
  
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
heatmap_func_ss <- function(Input_DF, COL_ANNOT, NoofRows, UNIT) {
  Heatmap_Data <- Input_DF
  Heatmap_Data <- Heatmap_Data[order(matrixStats:::rowVars(as.matrix(Heatmap_Data)), decreasing = TRUE), ]
  Heatmap_Data <- Heatmap_Data[c(1:NoofRows), ]
  Heatmap_Data <- Heatmap_Data[rowSums(is.na(Heatmap_Data)) != ncol(Heatmap_Data), ]
  Heatmap_Data <- as.data.frame(t(scale(t(Heatmap_Data))))
  Heatmap_Data <- as.data.frame(t(Heatmap_Data))
  
  # Sorting columns
  for(col in length(Heatmap_Data[1,]):1) {
    Column_Name = colnames(Heatmap_Data)[col]
    Heatmap_Data <- Heatmap_Data[order(Heatmap_Data[,Column_Name], decreasing = TRUE), ]
  }
  
  Heatmap_Data <- as.data.frame(t(Heatmap_Data))
  
  # Sorting rows
  for(col in length(Heatmap_Data[1,]):1) {
    Column_Name = colnames(Heatmap_Data)[col]
    Heatmap_Data <- Heatmap_Data[order(Heatmap_Data[,Column_Name], decreasing = TRUE), ]
  }
  
  REQD_SAMPLES <- as.character(intersect(rownames(COL_ANNOT), colnames(Heatmap_Data)))
  Heatmap_Annotation_Data <- COL_ANNOT[REQD_SAMPLES, ]
  Heatmap_Data <- Heatmap_Data[, REQD_SAMPLES]
  
  # Color sets
  n_sample_colors <- length(unique(Heatmap_Annotation_Data[, "Sample"]))
  n_tumor_colors <- length(unique(Heatmap_Annotation_Data[, "Tumor"]))
  
  Color_Set1 <- colorRampPalette(brewer.pal(min(9, n_sample_colors), "Set1"))(n_sample_colors)
  names(Color_Set1) <- unique(Heatmap_Annotation_Data[, "Sample"])
  
  Color_Set2 <- colorRampPalette(brewer.pal(min(8, n_tumor_colors), "Accent"))(n_tumor_colors)
  names(Color_Set2) <- unique(Heatmap_Annotation_Data[, "Tumor"])
  
  HEATMAP_ANNOTATION <- HeatmapAnnotation(
    Sample = as.character(Heatmap_Annotation_Data[, "Sample"]),
    Tumor = as.character(Heatmap_Annotation_Data[, "Tumor"]),
    col = list(Sample = Color_Set1, Tumor = Color_Set2)
  )
  
  col_fun = colorRamp2(c(min(Heatmap_Data, na.rm = TRUE), 0, max(Heatmap_Data, na.rm = TRUE)), c("blue", "white", "red"))
  
  Heatmap_Data <- Heatmap_Data[, rownames(Heatmap_Annotation_Data)]
  
  INPUT_DATA_HEATMAP <- Heatmap(
    as.matrix(Heatmap_Data),
    top_annotation = HEATMAP_ANNOTATION,    
    border_gp = gpar(col = "black", lwd = 2),
    name = as.character(UNIT),  # Ensure UNIT is a character
    column_split = Heatmap_Annotation_Data[, "Tumor"],
    rect_gp = gpar(col = "black", lwd = 0.2),
    row_gap = unit(5, "mm"),
    column_title = " ",
    column_title_side = "top",
    row_title = " ", 
    row_title_rot = 90,
    col = col_fun,
    row_names_gp = gpar(fontsize = 12, just = "center", fontface = "bold"),
    column_names_gp = gpar(fontsize = 14, just = "center", fontface = "bold"),
    column_title_gp = gpar(fontsize = 14, just = "center", fontface = "bold"),
    row_title_gp = gpar(fontsize = 14, just = "center", fontface = "bold"),
    column_names_side = "bottom",
    na_col = "grey90",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    column_names_centered = TRUE,
    column_names_max_height = unit(12, "cm"),
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_rot = 90,
    show_column_names = TRUE,
    heatmap_legend_param = list(direction = "vertical")
  )
  
  return(INPUT_DATA_HEATMAP)
}

######## Setup Project ----
## Initiate project
setupProject("Nanostring-LUCAT1") ; print(paste0("Working dir is: ", getwd()))
# If any project specific override: output folder if any
# output_dir <- paste0(output_dir, "/1.1_Nanostring")

## Install & Load Packages
cran_packages <- c("circlize", "clipr", "colorRamp2", "DT", "ggalluvial", "ggrepel", "grid", "igraph", "magick", "patchwork", "RColorBrewer", "tidyverse")
bioc_packages <- c("ComplexHeatmap","edgeR", "fgsea", "GSEABase", "GSVA", "limma", "msigdb", "msigdbr", "qusage", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
install_and_load_packages(cran_packages, bioc_packages)

######## Source & Process Input files ----
#Since the data is normalized counts; can not run the DESeq2 Pipeline
#So resorting to Limma Voom. Pipeline. Install edgeR, that installs limma as a dependency

# library("edgeR")
# library("tidyverse")
getwd()
file_samplesheet<-paste0(input_dir,"/samplesheet.csv")
file_source<-paste0(input_dir,"/LUCAT1 KD Initial Dataset .xlsx")
# file_gene_of_interest <- paste0(input_dir, "/Crispr screen gene targets.xlsx")

#Getting the Collected GeneSet of Interest from Excel List
# geneList <- readxl::read_excel(file_gene_of_interest, col_names = FALSE, sheet = "Sheet1")
# genesWithNa <- as.vector(as.matrix(geneList))
# genes <- values[!is.na(genesWithNa)] # Remove NA values
# GenesOfInterest <- unlist(genes) # Convert to simple list
# GenesOfInterest <- data.frame(GenesOfInterest)
# print(GenesOfInterest) # Print the list

#Sample Data for downstream analysis
Sample_Data<-read.csv(file = file_samplesheet, 
                      header=TRUE, 
                      stringsAsFactors = FALSE, #Made false to keep it as character
                      check.names = FALSE,
                      row.names = "TargetName")

# sampleAnnoFile from SegmentProperties sheet
# countFile & featureAnnoFile from BioProbeCountMatrix sheet
#### sampleAnnoFile----
#SegmentProperties Worksheet: SegmentDisplayName is default column
pre_sampleAnnoFile <- readxl::read_excel(file_source, 
                                         sheet="SegmentProperties")

next_sampleAnnoFile <- pre_sampleAnnoFile %>% 
  mutate(Samples = ifelse(Sample =="LUCAT1 KD", "LUCAT1_KD", Sample)) 
next_sampleAnnoFile <- next_sampleAnnoFile %>%
  mutate(Samples = ifelse(Samples =="CONTROL EDGE", "CONTROL", Samples))
next_sampleAnnoFile <- next_sampleAnnoFile %>%
  mutate(Samples = ifelse(Samples =="CONTROL INSIDE", "CONTROL", Samples))

next_sampleAnnoFile <- next_sampleAnnoFile %>% 
  mutate(Tumor = gsub("CONTROL EDGE", "TUMOR_EDGE", 
                      gsub("LUCAT1 KD", "TUMOR_WHOLE", 
                           gsub("CONTROL INSIDE", "TUMOR_INSIDE", Sample)))) 
sampleAnnoFile <- next_sampleAnnoFile %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#### featureAnnoFile----
#BioProbeCountMatrix Worksheet: TargetName is default column
pre_featureAnnoFile <- readxl::read_excel(file_source, 
                                          sheet = "BioProbeCountMatrix")

next_featureAnnoFile <- dplyr::select(pre_featureAnnoFile, 1:12) %>% #Get the featureData relevant columns only
  dplyr::select(TargetName, everything())

#Process the duplicates where the Targetname is not NegProbe-WTX and then combine them back
#Isolate rows with "NegProbe-WTX" in TargetName
featureAnnoFile_NegProbeWTX <- next_featureAnnoFile %>%
  filter(TargetName == "NegProbe-WTX")

featureAnnoFile_others <- next_featureAnnoFile %>%
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
  mutate(GeneID = next_featureAnnoFile$GeneID[match(TargetName, next_featureAnnoFile$TargetName)][1])

# Combine both processed datasets
featureAnnoFile <- bind_rows(featureAnnoFile_NegProbeWTX, featureAnnoFile_others) %>% 
  as.data.frame(., row.names = NULL, optional = FALSE, stringsAsFactors = FALSE)

#As TargetName (Which is the default Defined column of Nanostring Data) which holds gene Name and also multiple NegProbe-WTX 
#is not retained, so create anther column to retain it for DEG analysis.
featureAnnoFile[,"TargetGene"] <- featureAnnoFile[,"TargetName"]

#### countFile----
pre_countFile <- readxl::read_excel(file_source, 
                                    sheet="BioProbeCountMatrix") 
next_countFile <- pre_countFile %>%   dplyr::select(., c(3, 13:36)) %>% #Select the rest of sheet after the featureAnnoFile
  # filter(TargetName !=duplicates) %>%
  # setNames(c(gsub(" 1/2-", "", colnames(.)))) %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#Process the duplicates where the Targetname is not NegProbe-WTX and then combine them back
#countFile_NegProbeWTX <- countFile[which(countFile[,"TargetName"]=="NegProbe-WTX"),]
countFile_NegProbeWTX <- next_countFile %>% filter(TargetName == "NegProbe-WTX")

countFile_others <- next_countFile %>%
  filter(TargetName != "NegProbe-WTX") %>%
  # mutate_at(c("GeneID"), as.character) %>% 
  group_by(TargetName) %>%
  summarise(across(everything(), ~ sum(., na.rm = TRUE)), .groups = "drop")
# Combine
countFile <- bind_rows(countFile_NegProbeWTX, countFile_others) %>% 
  as.data.frame(., row.names = NULL, optional = FALSE, stringsAsFactors = FALSE)

######## CREATE SPATIAL EXPERIMENT OBJECT [seo]----
# BiocManager::install("standR")
# library(standR)
#install.packages("ggalluvial") #required for the standR package
#install.packages("magick") #required for the standR package

seo <- readGeoMx(countFile = countFile,
                 sampleAnnoFile = sampleAnnoFile,
                 featureAnnoFile = featureAnnoFile)

######## QC Steps ---- 
#Sample level QC
# library(ggplot2)
# library(ggalluvial)
# Create QC Directory
QC_DIR <- paste(output_dir, "/QC/", sep = "")
dir.create(QC_DIR)

#Visualize the data
graphics.off()
pdf(paste(QC_DIR, "Metadata.pdf", sep = ""), width=10, height = 10)
plotSampleInfo(seo, column2plot =c("SlideName", "Samples", "Tumor"))
graphics.off()

plotSampleInfo(seo, column2plot =c("Samples", "Tumor")) + coord_flip()

#### Gene level QC [seo_qc]----
seo #Dim 18676x24
names(colData(seo)) #46 Columns in ColData
view(seo@colData)
seo_qc <- addPerROIQC(seo, 
                      sample_fraction = 0.9, #Default
                      rm_genes =TRUE, #Default
                      min_count = 5) #Default
seo_qc #Gene with low count and expression values in more than threshold (sample_fraction=0.9)
# are removed by applying the function. Dim 19948x175, so removed 14 genes
dim(seo) 
dim(seo_qc) #Genes not meeting the above criteria were removed 19962 > 19948 

#If any genes are getting removed, we could plot those to visualize.
columns_of_interest <- c("SlideName")

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
plotGeneQC(seo_qc, top_n=12, ordannots = "Samples", col = Samples, point_size = 2)

#Numeric Values: AOISurfaceArea, AOINucleiCount, 
#RawReads, AlignedReads, DeduplicatedReads, TrimmedReads, StitchedReads, 
#SequencingSaturation, lib_size, countOfLowEprGenes, percentOfLowEprGene

#colData(seo)$regions
#colnames(seo) %>% print() 

#Final data of this segment: seo_qc
#### ROI level QC [seo_qc_roi]----
names(seo_qc@colData)
view(seo_qc@colData)
sapply(seo_qc@colData, class)
sapply(seo_qc@colData, is.numeric)

plotROIQC(seo_qc)

plotROIQC(seo_qc, 
          x_threshold = 25000,
          x_axis = "AOISurfaceArea",
          x_lab = "AOI Surface Area",
          # y_axis = "lib_size",
          # y_lab = "Library Size",
          y_axis = "RawReads",
          y_lab = "Raw Reads",
          y_threshold = 2.5e+06,
          color = Samples)

plotROIQC(seo_qc, 
          x_threshold = 25000, 
          x_axis = "AOISurfaceArea", 
          x_lab = "AOI Surface Area", 
          y_axis = "lib_size", 
          y_threshold = 2e+05,
          y_lab = "Library Size", 
          col = Samples)

plotROIQC(seo_qc, 
          x_threshold = 300, 
          x_axis = "AOINucleiCount", 
          x_lab = "AOI Nuclei Count", 
          y_axis = "RawReads", 
          # y_threshold = 2e+05,
          y_lab = "Raw Reads", 
          col = Samples)

colData(seo_qc)$AOINucleiCount 
#same as
seo_qc@colData$AOINucleiCount
view(seo_qc@colData)
#AOINuclei count of 150 looks like a good threshold from the figure
qc_keep <- colData(seo_qc)$AOINucleiCount > 300 &
  colData(seo_qc)$AOISurfaceArea > 25000 &
  colData(seo_qc)$RawReads > 2.5e+05 &
  colData(seo_qc)$AlignedReads > 2500000 &
  colData(seo_qc)$TrimmedReads > 80 &
  colData(seo_qc)$StitchedReads > 80 &
  colData(seo_qc)$SequencingSaturation > 50

table(qc_keep) # 3 Values Below threshold
dim(seo_qc) # Dim 18676x24
seo_qc_roi <- seo_qc[, qc_keep]
dim(seo_qc_roi) # Noothing Removed
# Comparing the Library Size with ROI Area size
view(seo_qc_roi@colData)
sum(seo_qc_roi@colData$Samples =="CONTROL")
sum(seo_qc_roi@colData$Samples =="LUCAT1_KD")

plotROIQC(seo_qc_roi, 
          x_threshold = 150, 
          x_axis = "AOINucleiCount", #Default
          y_threshold = 1e+01, 
          color = Samples)

plotROIQC(seo_qc_roi,
          x_axis = "SequencingSaturation",
          x_lab = "Sequencing Saturation",
          x_threshold = 50, 
          y_axis = "lib_size",
          y_lab = "Library Size",
          col = SlideName)


plotROIQC(seo_qc_roi,
          x_axis = "SequencingSaturation",
          x_lab = "Sequencing Saturation",
          y_axis = "AOISurfaceArea",
          y_lab = "AOISurfaceArea",
          col = SlideName)

plotROIQC(seo_qc_roi,
          x_axis = "RawReads",
          x_lab = "Raw Reads",
          x_threshold = 1e+06, 
          y_axis = "AOISurfaceArea",
          y_lab = "AOISurfaceArea",
          col = Samples)

plotROIQC(seo_qc_roi,
          x_axis = "RawReads",
          x_lab = "Raw Reads",
          x_threshold = 1e+06, 
          y_axis = "AOINucleiCount",
          y_lab = "AOINucleiCount",
          col = SlideName)

plotROIQC(seo_qc_roi,
          x_axis = "AlignedReads",
          x_lab = "Aligned Reads",
          x_threshold = 80,
          y_axis = "AOINucleiCount",
          y_lab = "AOINucleiCount",
          col = Tumor)

plotROIQC(seo_qc_roi,
          x_axis = "TrimmedReads",
          x_lab = "Trimmed Reads",
          x_threshold = 80,
          y_axis = "AOINucleiCount",
          y_lab = "AOINucleiCount",
          col = Tumor)

plotROIQC(seo_qc_roi,
          x_axis = "StitchedReads",
          x_lab = "Stitched Reads",
          x_threshold = 80,
          y_axis = "AOINucleiCount",
          y_lab = "AOINucleiCount",
          col = Samples)

# Relative log expression distribution
plotRLExpr(seo) #RLE of raw count 
plotRLExpr(seo_qc_roi)
#Remove the technical variations due to the library size differences
plotRLExpr(seo_qc_roi, ordannots = "Tumor", assay = 2, color = Tumor)+ ggtitle("Post QC Data")
#can also plot by tissue type or other classification
plotRLExpr(seo_qc_roi, ordannots = "Samples", assay = 2, color = Samples)+ ggtitle("Post QC Data")


######## Dimentionality Reduction [seoPCA->seoUMAP]----
#seo_qc_roi@assays #Assay 2 is based on logcounts
# BiocManager::install("scater")
drawPCA(seo_qc_roi, assay = 2, color = Samples) # however since the pca will change axis every time we plot, 
#We can save the data to analyze it the same way every time
drawPCA(seo_qc_roi, assay =2, color = Tumor)
#To make it reproducible
set.seed(100)

seoPCA <-  scater::runPCA(seo_qc_roi)
#runPCA adds reducedDimNames PCA
pca_results <-  reducedDim(seoPCA, "PCA")
drawPCA(seoPCA, precomputed = pca_results, col = Tumor)
drawPCA(seoPCA, precomputed = pca_results, col = Samples)

#Draw PCA Scree Plot
plotScreePCA(seo_qc_roi, precomputed = pca_results)
#Plot Pair PCA
plotPairPCA(seo_qc_roi, col= tissue, precomputed = pca_results, n_dimension = 4)
plotPairPCA(seo_qc_roi, col= SlideName, precomputed = pca_results, n_dimension = 4)

#Plot Multidimensional Scaling (MDS) plot
standR::plotMDS(seo_qc_roi, assay = 2, color = Samples)
standR::plotMDS(seo_qc_roi, assay = 2, color = Tumor)


#UMAP
set.seed(100)

seoUMAP <- scater::runUMAP(seoPCA, dimred = "PCA")
#runUMAP adds one more reducedDimNames UMAP
plotDR(seoUMAP, dimred = "UMAP", col = Samples)
plotDR(seoUMAP, dimred = "UMAP", col = Tumor)

names(seoUMAP@metadata)

######## ↓ ↓ Q3 Normalization and Saving Data [seoUMAP -> seoUMAP_Normalized_Q3] ----
# Input: seoUMAP; output: seoUMAP_Normalized_Q3
#Make Count Data from Source seoUMAP----
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
##(SKIP ENDS)
#### Add Normalized Data to logcounts field----
NORMALIZATION <- "upperquartile"

seoUMAP_Normalized_Q3 <- geomxNorm(seoUMAP, method=NORMALIZATION, log = T)

names(seoUMAP_Normalized_Q3@metadata) #Added "norm.method" and "norm.factor" to metadata
seoUMAP_Normalized_Q3@metadata$norm.method
seoUMAP_Normalized_Q3@metadata$norm.factor

plotRLExpr(seoUMAP_Normalized_Q3, assay = 2, color = Samples) + ggtitle("Q3 Normalized")
plotRLExpr(seoUMAP_Normalized_Q3, assay = 2, color = Tumor) + ggtitle("Q3 Normalized")
##Make Normalized Count Data ----

Normalized_Counts_Data_Q3 <- as.data.frame(logcounts(seoUMAP_Normalized_Q3))

##(SKIP IF NOT SAVING Normalized Count Data)Generate a Gene Column for Data Saving ----
Normalized_Counts_Data_Q3 <- as.data.frame(cbind(Gene=rownames(Normalized_Counts_Data), Normalized_Counts_Data))
write.table(Normalized_Counts_Data_Q3, paste(output_dir, NORMALIZATION, "_Normalized_Counts.txt", sep=""), 
            sep="\t", row.names = F, col.names = T, quote = T)
#Read Data
Normalized_Counts_Data_Q3 <- read.delim(paste(output_dir, NORMALIZATION, "_Normalized_Counts.txt", sep=""), 
                                        sep="\t", header = T, check.names = F, stringsAsFactors = F)
#Assign the rownames back as Genes
rownames(Normalized_Counts_Data_Q3) <- Normalized_Counts_Data_Q3[,"Gene"]

#Remove the Gene Column to make original normalized counts
Normalized_Counts_Data_Q3[,"Gene"] <- NULL
##(SKIP ENDS)
##Make Feature Data ----
Feature_Data <- as.data.frame(rowData(seoUMAP_Normalized_Q3))
#Skip if not saving
# Feature_Data <- as.data.frame(cbind(Gene=rownames(Feature_Data), Feature_Data))
##Make Sample Data ----
Sample_Data <- as.data.frame(colData(seoUMAP_Normalized_Q3))
#Skip if not saving
#Sample_Data <- as.data.frame(cbind(ROI=rownames(Sample_Data), Sample_Data))


#↑ ↑ Q3 Normalization End~~~~~~~~~~~~~~~~~~

######## ↓ ↓ TMM Normalization and Saving Data [seoUMAP -> seoUMAP_Normalized_TMM]----
# Input: seoUMAP output: seoUMAP_Normalized_TMM
#Data has become seoUMAP now with addition of two reducedDimNames PCA and UMAP
#names(seoUMAP@metadata)
NORMALIZATION_METHOD <- "TMM"
seoUMAP_Normalized_TMM <- geomxNorm(seoUMAP, method = NORMALIZATION_METHOD, log = TRUE) #TMM Normalization
#geomxNorm Adds two more items to metadata: norm.factor, norm.method

#names(seoUMAP_Normalized_TMM@metadata)
plotRLExpr(seoUMAP_Normalized_TMM, assay = 2, color = Samples) + ggtitle("TMM Normalized")
plotRLExpr(seoUMAP_Normalized_TMM, assay = 2, color = Tumor) + ggtitle("TMM Normalized")

seoUMAP_Normalized_TMM@metadata$norm.method
seoUMAP_Normalized_TMM@metadata$norm.factor
#Make Count Data from Source seoUMAP----
Counts_Data <- as.data.frame(counts(seoUMAP))
#Make Normalized Count from TMM----
Normalized_Counts_Data_TMM <- as.data.frame(logcounts(seoUMAP_Normalized_TMM))

###Make Sample Data -Same from Q3 or TMM File----
Sample_Data <- as.data.frame(colData(seoUMAP_Normalized_TMM))
##Make Feature Data ----
Feature_Data <- as.data.frame(rowData(seoUMAP_Normalized_TMM))
##↓↓~~~~~~~~~~~~~~~~~~~~~

####LOOP Plot QC Figures in Loop----
names(seo_qc_roi@colData)
#Inputs seoUMAP_Normalized_TMM, seoUMAP_Normalized_Q3 & gmt files in input/MSIG_DB
ROI_ANNOTATION_COLS <- c("Sample", "Samples", "Tumor")
dictionary <- c(TMM=seoUMAP_Normalized_TMM, Q3=seoUMAP_Normalized_Q3) #seoUMAP_Normalized_TMM (for TMM) or seoUMAP_Normalized_Q3 (for Q3)
# eval(GRP)
i = 1
NORM = "Q3" #[INPUT_NEEDED] TMM or Q3 for naming the plots and files

for(i in 1:length(ROI_ANNOTATION_COLS)){
  GRP <- ROI_ANNOTATION_COLS[[i]]
  print(GRP)
  set.seed(100)
  RLE_PLOT <- plotRLExpr(dictionary[[NORM]], assay = 2, ordannots=GRP, color=get(GRP)) + ggtitle(paste(NORM, GRP, "RLE")) + labs(color=GRP)
  PAIRWISE_PCA <- plotPairPCA(dictionary[[NORM]], title=GRP, col = get(GRP), shape = get(GRP), assay=2, n_dimension=2) + labs(color=GRP, shape=GRP)
  PCA_BI_PLOT <- plotPCAbiplot(dictionary[[NORM]], n_loadings=10, assay=2, col=get(GRP)) + labs(color=GRP)
  MultiDIM_SCALING_PLOT <- standR::plotMDS(dictionary[[NORM]], assay=2, col=get(GRP), shape=get(GRP)) + labs(color=GRP, shape=GRP)
  UMAP <- plotDR(dictionary[[NORM]], dimred="UMAP", col=get(GRP)) + labs(color=GRP)
  PCA_SCREE_PLOT <- plotScreePCA(dictionary[[NORM]], assay=2, dims=10)
  
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
# library("ComplexHeatmap")
# library("GSVA")
# library("RColorBrewer")
# install.packages("colorRamp2")
# library("colorRamp2")
# BiocManager::install("qusage")
# library("qusage")
Normalized_Counts_Data <- Normalized_Counts_Data_TMM

MSIG_DB <- paste0(input_dir, "/MSIG_DB/")
GeneSets <- list.files(MSIG_DB, pattern= "*.gmt", full.names = F)
MSigDB_Dictionary <- list(h_all="HallMark Gene Sets", 
                          c2_cp_biocarta="BioCarta subset of Canonical Pathways",
                          c2_cp_pid="PID subset of Canonical Pathways",
                          c2_cp_reactome="Reactome subset of Canonical Pathways",
                          c2_cp_wikipathways="WikiPathways subset of Canonical Pathways")

geneset =1
for (geneset in 1:length(GeneSets)){
  # geneset_name <-  "c2_cp_biocarta" #For Non-loop version
  geneset_name = gsub("\\.","_",
                      gsub(".v2023.2.Hs.symbols.gmt", "", GeneSets[geneset]))
  
  print(geneset_name)
  ####
  gset_human=paste(MSIG_DB, GeneSets[geneset], sep = "")
  # Signature <- qusage::read.gmt(paste(MSIG_DB,"c2.cp.biocarta.v2023.2.Hs.symbols.gmt", sep="")) #Non-loop version
  Signature <- qusage::read.gmt(gset_human)
  ####
  
  ssParam <- gsvaParam(as.matrix(Normalized_Counts_Data),
                       Signature,
                       kcdf = "Gaussian",
                       maxDiff = TRUE,
                       minSize = 2)  # *NEW from last analysis: Set minimum gene set size to 2
  
  ssGSEAScores <- GSVA::gsva(ssParam, verbose=TRUE)
  
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
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("PID_","",as.character(x)))
  rownames(Input_DF) <- sapply(rownames(Input_DF),function(x) gsub("WP_","",as.character(x)))
  ###########################
  # COL_ANNOT <- NULL
  # COL_ANNOT = Sample_Data[,c("ROI","ScanLabel","Neuron","tissue","CD45")]
  # COL_ANNOT[,"Group"] <- COL_ANNOT[,"ScanLabel"]
  # rownames(COL_ANNOT) <- COL_ANNOT[,"ROI"]
  # COL_ANNOT[,"ROI"] <- NULL
  COL_ANNOT <- Sample_Data
  ##########################
  No_Of_Rows = 100
  # UNIT = geneset_name
  UNIT = get("MSigDB_Dictionary")[[geneset_name]]
  #########################3
  
  ssGSEA_Heatmap <- heatmap_func_ss(Input_DF,COL_ANNOT,100, UNIT) #*NEW from last analsysis, added UNIT, improved the function
  
  #NCOLS = as.character(length(Heatmap_Data[,1]))
  #NROWS = as.character(length(Heatmap_Data[1,]))
  #######################################
  graphics.off()
  # pdf_height <- max(30, nrow(Input_DF) * 0.2)  # Adjust the multiplier as needed
  pdf_height <- 30
  pdf(paste(output_dir, "/", geneset_name,"_ssGSEAScores.pdf",sep=""),width =60,height =pdf_height) #Changed to variable instead of 30
  draw(ssGSEA_Heatmap,padding = unit(c(1,5,1,1), "in"),heatmap_legend_side = "right",row_title = "", row_title_gp = gpar(col = "red"),legend_grouping = "original",
       column_title = paste("ssGSEA- ", UNIT,sep=""), column_title_gp = gpar(fontsize = 32))
  graphics.off()
  #####################################
}
