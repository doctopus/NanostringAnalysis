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
  
  Color_Set1 <- Color_Sets[[1]][c(1:length(as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Sample"]),"Sample"])))]
  names(Color_Set1) <- as.character(Heatmap_Annotation_Data[!duplicated(Heatmap_Annotation_Data[,"Sample"]),"Sample"])
  
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
  
  
  
  HEATMAP_ANNOTAION <- HeatmapAnnotation(Sample=as.character(Heatmap_Annotation_Data[,"Sample"]),
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
setupProject("Nanostring-ALPK3") ; print(paste0("Working dir is: ", getwd()))
# If any project specific override: output folder if any
# output_dir <- paste0(output_dir, "/1.1_Nanostring")

## Install & Load Packages
cran_packages <- c("colorRamp2", "DT", "ggalluvial", "ggrepel", "igraph", "magick", "RColorBrewer", "tidyverse")
bioc_packages <- c("ComplexHeatmap", "edgeR", "GSEABase", "GSVA", "limma", "msigdb", "qusage", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
install_and_load_packages(cran_packages, bioc_packages)

######## Source & Process Input files ----
# sampleAnnoFile from SegmentProperties sheet
# countFile & featureAnnoFile from BioProbeCountMatrix sheet
file_samplesheet<-paste0(input_dir,"/samplesheet_2.csv")
file_source<-paste0(input_dir,"/ALPK3 experiment Initial Dataset .xlsx")
file_gene_of_interest <- paste0(input_dir, "/Crispr screen gene targets.xlsx")
#### sampleAnnoFile----
#SegmentProperties Worksheet: SegmentDisplayName is default column
pre_sampleAnnoFile <- readxl::read_excel(file_source, sheet="SegmentProperties")

# colnames(pre_sampleAnnoFile)[colnames(pre_sampleAnnoFile)=='sample'] <- 'Sample'
# don't rename column "sample" as we'll remove it altogether after. getting data from it
sampleAnnoFile <- pre_sampleAnnoFile %>%
  mutate(
    Sample = case_when(
      is.na(sample) ~ NA_character_,
      grepl("ALPK3 KD", sample) ~ "ALPK3 KD",
      grepl("Control", sample) ~ "Control",
      TRUE ~ NA_character_
    ),
    Tumor = case_when(
      is.na(sample) ~ NA_character_,
      grepl("Tumor Edge", sample) ~ "Tumor Edge",
      grepl("Tumor Inside", sample) ~ "Tumor Inside",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-sample) %>% 
  as.data.frame(., row.names=NULL, optional=FALSE, stringAsFactors = FALSE)

#### featureAnnoFile----
#BioProbeCountMatrix Worksheet: TargetName is default column
pre_featureAnnoFile <- readxl::read_excel(file_source, sheet = "BioProbeCountMatrix")

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
pre_countFile <- readxl::read_excel(file_source,sheet="BioProbeCountMatrix") 
next_countFile <- pre_countFile %>%   dplyr::select(., c(3, 13:79)) %>% #Select the rest of sheet after the featureAnnoFile
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

#SampleData in-case Needed in Downstream----
Sample_Data<-read.csv(file = file_samplesheet,
                     header=TRUE,
                     stringsAsFactors = TRUE, #Make false to keep it as character
                     check.names = FALSE,
                     row.names = "TargetName")

#Getting the Collected GeneSet of Interest from Excel List for Downstream ----
geneList <- readxl::read_excel(file_gene_of_interest, col_names = FALSE, sheet = "Sheet1")
genesWithNa <- as.vector(as.matrix(geneList))
genes <- values[!is.na(genesWithNa)] # Remove NA values
GenesOfInterest <- unlist(genes) # Convert to simple list
GenesOfInterest <- data.frame(GenesOfInterest)
# print(GenesOfInterest) # Print the list

######## CREATE SPATIAL EXPERIMENT OBJECT [seo]----
seo <- readGeoMx(countFile = countFile,
                 sampleAnnoFile = sampleAnnoFile,
                 featureAnnoFile = featureAnnoFile)

######## QC Steps ---- 
#Sample level QC
# library(ggplot2)
# library(ggalluvial)
#Visualize the data
plotSampleInfo(seo, column2plot =c("Sample", "Tumor"))
plotSampleInfo(seo, column2plot =c("SlideName", "Sample", "Tumor")) + coord_flip()
#### Gene level QC [seo_qc]----
seo #Dim 18676x67
names(colData(seo)) #30 Columns in ColData
seo_qc <- addPerROIQC(seo, 
                      sample_fraction = 0.9, #Default
                      rm_genes =TRUE, #Default
                      min_count = 5) #Default
seo_qc #Gene with low count and expression values in more than threshold (sample_fraction=0.9)
# are removed by applying the function. Dim 19948x175, so removed 14 genes
dim(seo) 
dim(seo_qc) #Genes not meeting the above criteria were removed [Nothing] 18676x67 > 18676x67
names(seo_qc@colData)
seo_qc@colData
view(colData(seo_qc)[ , c("countOfLowEprGene", "percentOfLowEprGene", "ScanLabel", "lib_size", "tissue")])
view(colData(seo_qc))

length(seo@metadata$genes_rm_rawCount) #0
length(seo_qc@metadata$genes_rm_rawCount) #67 Data added. to the field
names(seo_qc@metadata)
metadata(seo_qc) |> names() #Same as above
metadata(seo) |>names()

#addPerROIQC added columns to colData :lib_size, countOfLowEprGene, percentOfLowEprGene
# and also added columns to metadata: lcpm_threshold, genes_rm_rawCount, genes_rm_logCPM  


# plotGeneQC(seo_qc, ordannots = "regions", col = regions, point_size = 2)
plotGeneQC(seo_qc)

plotGeneQC(seo_qc, top_n=12, ordannots = "SlideName", col = SlideName, point_size = 2)
plotGeneQC(seo_qc, top_n=12, ordannots = "Sample", col = Sample, point_size = 2)
plotGeneQC(seo_qc, top_n =12, ordannots = "NF-H-_CD45-_Edge_Gland", col = 'NF-H-_CD45-_Edge_Gland', point_size = 2 )

sapply(seo_qc@colData, class)
#Numeric Values: AOISurfaceArea, AOINucleiCount, 
#RawReads, AlignedReads, DeduplicatedReads, TrimmedReads, StitchedReads, 
#SequencingSaturation, lib_size, countOfLowEprGenes, percentOfLowEprGene

#colData(seo)$regions
#colnames(seo) %>% print() 

#Final data of this segment: seo_qc
#### ROI level QC [seo_qc_roi]----
names(colData(seo_qc))

plotROIQC(seo_qc)

plotROIQC(seo_qc, x_threshold = 900, color = Sample)
plotROIQC(seo_qc, 
          x_axis = "lib_size", x_threshold = 1e+01, 
          y_axis = "AOISurfaceArea", y_threshold = "AOISurfaceArea",
          col = Sample)

colData(seo_qc)$AOINucleiCount 
#same as
seo_qc@colData$AOINucleiCount

#AOINuclei count of 150 looks like a good threshold from the figure
qc <- colData(seo_qc)$AOINucleiCount > 900 
table(qc) # 3 Values Below threshold
dim(seo_qc) # Dim 19948x175
seo_qc_roi <- seo_qc[, qc]
dim(seo_qc_roi) # We removed 3 ROI (samples/columns)from dataset. Dim 18676 >
# Comparing the Library Size with ROI Area size

plotROIQC(seo_qc, 
          x_threshold = 20000, 
          x_axis = "AOISurfaceArea", 
          x_lab = "AreaSize", 
          y_axis = "lib_size", 
          y_lab = "Library Size", 
          col = Sample)

plotROIQC(seo_qc_roi, 
          # x_threshold = 150, 
          x_axis = "AOINucleiCount", #Default
          # y_threshold = 1e+01, 
          color = Sample)

plotROIQC(seo_qc_roi,
          x_axis = "SequencingSaturation",
          x_lab = "Sequencing Saturation",
          y_axis = "lib_size",
          y_lab = "Library Size",
          col = SlideName)


plotROIQC(seo_qc_roi,
          x_axis = "SequencingSaturation",
          x_lab = "Sequencing Saturation",
          y_axis = "AOISurfaceArea",
          y_lab = "AOISurfaceArea",
          col = SlideName)

# Relative log expression distribution
plotRLExpr(seo) #RLE of raw count 
plotRLExpr(seo_qc_roi)
#Remove the technical variations due to the library size differences
plotRLExpr(seo_qc_roi, ordannots = "SlideName", assay = 2, color = SlideName)
#can also plot by tissue type or other classification
plotRLExpr(seo_qc_roi, ordannots = "tissue", assay = 2, color = tissue)

#### ssGSEA Analysis----
# library("ComplexHeatmap")
# library("GSVA")
# library("RColorBrewer")
# install.packages("colorRamp2")
# library("colorRamp2")
# BiocManager::install("qusage")
# library("qusage")

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
  UNIT = get("MSigDB_Dictionary")[geneset_name]
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

##↑ ↑~~~~~~~~~~~~~~~~~~~~~~ ----

COMPARISIONS <- as.data.frame(cbind(COMPARE_GROUP_NAME=c("ScanLabel","ScanLabel","ScanLabel_Neuron","ScanLabel_Neuron","ScanLabel_Neuron","ScanLabel_Neuron"),
                                    GROUP1_NAME=c("Sympa","FP","Sympa_NF_H_POS","FP_NF_H_POS","Sympa_NF_H_POS","Sympa_NF_H_NEG"),
                                    GROUP2_NAME=c("FP","BBP","Sympa_NF_H_NEG","FP_NF_H_NEG","FP_NF_H_POS","FP_NF_H_NEG")))
#
comp = 1
for(comp in 1:length(COMPARISIONS[,1])){
  print(comp)
  COMPARE_GROUP_NAME = COMPARISIONS[comp,"COMPARE_GROUP_NAME"] 
  GROUP1_NAME = COMPARISIONS[comp,"GROUP1_NAME"] 
  GROUP2_NAME <-  COMPARISIONS[comp,"GROUP2_NAME"]
  
  ###################################################################  
  print(paste(COMPARE_GROUP_NAME,":",GROUP1_NAME,"_Vs_",GROUP2_NAME))
  Comparision_Sample_Data <- Experiment_Sample_Data[,c("ROI","ScanLabel","Neuron","tissue","CD45")]
  Comparision_Sample_Data[,"Group"] <- Experiment_Sample_Data[,COMPARE_GROUP_NAME]
  print(table(Comparision_Sample_Data[,"Group"]))
}

#-----
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


#countFile----
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
