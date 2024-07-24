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

##### Setup Project ----
## Initiate project
setupProject("Nanostring-ALPK3") ; print(paste0("Working dir is: ", getwd()))
# If any project specific override: output folder if any
# output_dir <- paste0(output_dir, "/1.1_Nanostring")

## Install & Load Packages
cran_packages <- c("colorRamp2", "DT", "ggalluvial", "ggrepel", "igraph", "magick", "RColorBrewer", "tidyverse")
bioc_packages <- c("ComplexHeatmap", "edgeR", "GSEABase", "GSVA", "limma", "msigdb", "qusage", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
install_and_load_packages(cran_packages, bioc_packages)

#### Source & Process Input files ----
#Since the data is normalized counts; can not run the DESeq2 Pipeline
#So resorting to Limma Voom. Pipeline. Install edgeR, that installs limma as a dependency

# library("edgeR")
# library("tidyverse")
getwd()
file_samplesheet<-paste0(input_dir,"/samplesheet_2.csv")
file_source<-paste0(input_dir,"/ALPK3 KD in 358 Nanostring Q3 Normalized Data.xlsx")
file_gene_of_interest <- paste0(input_dir, "/Crispr screen gene targets.xlsx")

## Process input files

normalized_counts <- readxl::read_excel(file_source, sheet = "TargetCountMatrix")
Normalized_Counts_Data <- data.frame(normalized_counts, check.names= FALSE, row.names = "TargetName")

featureData <- readxl::read_excel(file_source, col_types= "text", sheet = "BioProbeProperties")
featureData <- featureData %>% filter(TargetName!="NegProbe-WTX")
# ↓ Does not work since NegProbe-WTX are many. In normalized counts however, it is one. 
# May need to remove or deal with these NegProbe-WTX to use Features_Data 
# Removed all where TargetName is NegProbe-WTX
Features_Data <- data.frame(featureData, check.names = FALSE, row.names = "TargetName")

Sample_Data<-read.csv(file = file_samplesheet, 
                     header=TRUE, 
                     stringsAsFactors = FALSE, #Made false to keep it as character
                     check.names = FALSE,
                     row.names = "TargetName")

#Getting the Collected GeneSet of Interest from Excel List
geneList <- readxl::read_excel(file_gene_of_interest, col_names = FALSE, sheet = "Sheet1")
genesWithNa <- as.vector(as.matrix(geneList))
genes <- values[!is.na(genesWithNa)] # Remove NA values
GenesOfInterest <- unlist(genes) # Convert to simple list
GenesOfInterest <- data.frame(GenesOfInterest)
# print(GenesOfInterest) # Print the list

#Old Code for reference (May Delete)----
# sampleData<-read.csv(file = file_samplesheet, 
#                      header=TRUE, 
#                      stringsAsFactors = FALSE, #Need as characters for this heatmap code function
#                      check.names = FALSE,
#                      row.names = "sample")
# sampleData <- sampleData[,c("cellLine", "drug")]
# sampleData <- sampleData %>% mutate(drug = gsub("-", ".", drug))
# sampleData[, c("cellLine", "drug")] <- lapply(sampleData[, c("cellLine", "drug")], factor)
# head(sampleData)

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
