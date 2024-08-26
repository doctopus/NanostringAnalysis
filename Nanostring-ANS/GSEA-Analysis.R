#GSVA Analysis and Heatmap creation
#For Custom Gene Sets such as WNT, NOTCH, HEDGEHOG pathway related genes.
#Pathways to focus on could be filtered from specific genesets such as Hallmark, Reactome, Biocarta etc
#Input: Normalized Count Data & Sample Data
#Output: Heatmap Plots
library(msigdbr)

#May SKIP Setting Up Environment if Done already-----
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
## Install & Load Packages
cran_packages <- c("circlize", "clipr", "colorRamp2", "DT", "ggalluvial", "ggrepel", "grid", "igraph", "magick", "patchwork", "RColorBrewer", "tidyverse")
bioc_packages <- c("ComplexHeatmap","edgeR", "fgsea", "GSEABase", "GSVA", "limma", "msigdb", "msigdbr", "qusage", "SpatialExperiment", "SpatialDecon", "speckle", "standR", "vissE")
install_and_load_packages(cran_packages, bioc_packages)

#### Non-Function Way of Getting Pathway specific genes from MSIGDB----
M <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
M2 <- msigdbr(species = "Mus musculus", category = "C2")
#Unique number of genes M
print(length(unique(M$gene_symbol))) #1274
print(length(unique(M$gs_name))) #29

#Unique number of genes M2
print(length(unique(M2$gene_symbol))) #17616
print(length(unique(M2$gs_name))) #6366

####. Finding Genes of a category of. pathway from msigdb ---- 
Human_Hallmark <- msigdbr(species="Homo sapiens", category ="H")
MTORC1_Signaling_Pathway_Genes <- Human_Hallmark %>% dplyr::filter(gs_name =="HALLMARK_MTORC1_SIGNALING") %>% 
  dplyr::select(c("gs_name", "gene_symbol", "gs_description"))
# install.packages("clipr")
# library(clipr)
# Export the data frame to clipboard
write_clip(MTORC1_Signaling_Pathway_Genes)


#### Retrieve mouse gene sets ----
msigdb_mouse = msigdbr(species = "Mus musculus") #3783805

#Unique number of genes
print(length(unique(msigdb_mouse$gene_symbol))) #17961
print(length(unique(msigdb_mouse$gs_name))) #32872

# Filter for Wnt signaling-related gene sets
wnt_genesets_mouse = msigdb_mouse %>%
  filter(grepl("WNT", gs_name, ignore.case = TRUE) | 
           grepl("WNT", gs_description, ignore.case = TRUE))

#Filter Immune related pathways in All Mouse Genesets
immune_pathways_mouse = msigdb_mouse %>% 
  filter(grepl("IMMUN", gs_name, ignore.case = TRUE) |
           grepl("IMMUN", gs_description, ignore.case = TRUE))
print(length(unique(immune_pathways_mouse$gs_name))) #346
print(length(unique(immune_pathways_mouse$gene_symbol))) #9427

#immune_pathways_unique <- unique(immune_pathways_mouse$gs_name)
immune_pathways_genes_unique <- unique(immune_pathways_mouse$gene_symbol)

head(immune_pathways_genes_unique)

#Mouse Immune Canonical C2
immune_pathways_mouse_c2 = M2 %>% 
  filter(grepl("IMMUN", gs_name, ignore.case = TRUE) |
           grepl("IMMUN", gs_description, ignore.case = TRUE))
print(length(unique(immune_pathways_mouse_c2$gs_name))) #54
print(length(unique(immune_pathways_mouse_c2$gene_symbol))) #3341
immune_pathways_genes_unique_c2 <- unique(immune_pathways_mouse_c2$gene_symbol)

#Mouse Immune Canonical C2:CP
immune_pathways_mouse_c2cp = M %>% 
  filter(grepl("IMMUN", gs_name, ignore.case = TRUE) |
           grepl("IMMUN", gs_description, ignore.case = TRUE))
print(length(unique(immune_pathways_mouse_c2cp$gs_name))) #0
print(length(unique(immune_pathways_mouse_c2cp$gene_symbol))) #0
immune_pathways_genes_unique_c2cp <- unique(immune_pathways_mouse_c2cp$gene_symbol)



immune_pathways_df <- data.frame("Immune Related Pathways" = unlist(immune_pathways_unique))
write_clip(immune_pathways_df)
# Display the data frame
print(immune_pathways_df)


#### Extract unique genes associated with Wnt signaling ----
wnt_genes_mouse = unique(wnt_genesets_mouse$gene_symbol)
#Check if a particular gene of interest is in the list
"Dkk4" %in% wnt_genes_mouse
# Print the number of Wnt-related gene sets and genes
cat("Number of Wnt-related gene sets:", nrow(wnt_genesets_mouse), "\n") #8058
cat("Number of unique Wnt-related genes:", length(wnt_genes_mouse), "\n") #3724

# To see the gene set names:
print(length(unique(wnt_genesets_mouse$gs_name)))

# To focus on a specific gene set, e.g., REACTOME_SIGNALING_BY_WNT:
reactome_wnt_genes_mouse = wnt_genesets_mouse %>%
  filter(gs_name == "REACTOME_SIGNALING_BY_WNT") %>%
  pull(gene_symbol) %>%
  unique()

cat("Number of genes in REACTOME_SIGNALING_BY_WNT:", length(reactome_wnt_genes_mouse), "\n") #316

#### Na Channel Reelated Genes ----
msigdb_mouse = msigdbr(species = "Mus musculus") #3783805

#Unique number of genes
print(length(unique(msigdb_mouse$gene_symbol))) #17961
print(length(unique(msigdb_mouse$gs_name))) #32872

#Filter Na Channel related pathways in All Mouse Genesets
naChannel_pathways_mouse = msigdb_mouse %>% 
  filter((grepl("SODIUM", gs_name, ignore.case = TRUE) |
           grepl("SODIUM", gs_description, ignore.case = TRUE)) &
         (grepl("CHANNEL", gs_name, ignore.case = TRUE) |
           grepl("CHANNEL", gs_description, ignore.case = TRUE)))
#Same result as below
naChannel_pathways_mouse = msigdb_mouse %>% 
  filter((grepl("SODIUM", gs_description, ignore.case = TRUE)) &
          (grepl("CHANNEL", gs_description, ignore.case = TRUE)))

print(length(unique(naChannel_pathways_mouse$gs_name))) #13
print(length(unique(naChannel_pathways_mouse$gene_symbol))) #102

naChannel_pathways_unique <- unique(naChannel_pathways_mouse$gs_name)
naChannel_pathways_genes_unique <- unique(naChannel_pathways_mouse$gene_symbol)

head(naChannel_pathways_genes_unique)
summary(naChannel_pathways_genes_unique)

head(naChannel_pathways_unique)
summary(naChannel_pathways_unique)

#### Useful Functions ----

## To Find genes of a particular pathway from its name
find_genes_for_specific_pathway <- function(pathway_name, species = "Mus musculus") {
  # Get all gene sets for the species
  all_gene_sets <- msigdbr(species = species)
  
  # Filter for the specific pathway
  pathway_genes <- all_gene_sets %>%
    filter(gs_name == pathway_name) %>%
    pull(gene_symbol) %>%
    unique()
  
  return(pathway_genes)
}
# Example usage:
pathway_of_interest <- "REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING"
species <- "Mus musculus"
genes_in_pathway <- find_genes_for_specific_pathway(pathway_of_interest, species)
print(paste("Number of genes in", pathway_of_interest, ":", length(genes_in_pathway)))
print("Genes:")
print(genes_in_pathway)


## To check which collection the particular pathway belongs to as well; it gives both the genes and the collection, eg C2 that the pathway belongs to.
find_genes_and_collection_for_specific_pathway <- function(pathway_name, species = "Mus musculus") {
  # Get all gene sets for the species
  all_gene_sets <- msigdbr(species = species)
  # Filter for the specific pathway
  pathway_info <- all_gene_sets %>%
    filter(gs_name == pathway_name)
  
  if(nrow(pathway_info) == 0) {
    return(list(genes = character(0), collection = NA))
  }
  
  genes <- unique(pathway_info$gene_symbol)
  collection <- unique(pathway_info$gs_cat)
  return(list(genes = genes, collection = collection))
}
# Example usage:
pathway_of_interest <- "REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING"
species <- "Mus musculus"
result <- find_genes_and_collection_for_specific_pathway(pathway_of_interest, species)
print(paste("Pathway:", pathway_of_interest))
print(paste("Collection:", result$collection))
print(paste("Number of genes:", length(result$genes)))
print("Genes:")
print(result$genes)












#### Functions Way of Code ----

# Function to get gene sets for a specific pathway and collection
get_pathway_genesets <- function(species, pathway, collections) {
  msigdb_data <- msigdbr(species = species)
  
  pathway_genesets <- msigdb_data %>%
    filter(
      (gs_cat %in% c("H", "C2")) &
        (gs_subcat %in% c("", "CP:REACTOME", "CP:BIOCARTA", "CP:WIKIPATHWAYS")) &
        (grepl(pathway, gs_name, ignore.case = TRUE) | 
           grepl(pathway, gs_description, ignore.case = TRUE)) &
        (gs_cat == "H" | gs_subcat %in% paste0("CP:", toupper(collections)))
    ) %>%
    split(.$gs_name) %>%
    lapply(function(x) unique(x$gene_symbol))
  
  return(pathway_genesets)
}

# Function to perform GSVA and return enrichment scores
perform_gsva <- function(expr_data, gene_sets) {
  gsvaParam <- gsvaParam(exprData = as.matrix(expr_data),
                         geneSets = gene_sets,
                         kcdf = "Gaussian")
  gsva_result <- gsva(gsvaParam, verbose = TRUE)
  return(gsva_result)
}

# Main analysis function
analyze_pathway <- function(expr_data, sample_data, species, pathway, collections) {
  # Get gene sets
  gene_sets <- get_pathway_genesets(species, pathway, collections)
  
  # Perform GSVA
  gsva_scores <- perform_gsva(expr_data, gene_sets)
  
  # Return only the GSVA scores
  return(gsva_scores)
}


## Load data----
Normalized_Counts_Data <- Normalized_Counts_Data_Q3 #Or Normalized_Counts_Data_TMM
Sample_Data <- Sample_Data %>% mutate(SlideName = gsub("BBP", "BFP", SlideName))
# Ensure that the column names of Normalized_Counts_Data match the row names of Sample_Data
if (all(colnames(Normalized_Counts_Data_TMM) == rownames(Sample_Data))) {
  message("Normalized Data and Sample Data Match")
} else {
  stop("Sample names in Normalized_Counts_Data and Sample_Data do not match")
}

# Define variables
species <- "Mus musculus"
pathways <- c("WNT", "HEDGEHOG", "NOTCH")
collections <- c("HALLMARK", "BIOCARTA", "REACTOME", "WIKIPATHWAYS")

# Perform analysis for each pathway
results_list <- lapply(pathways, function(pathway) {
  analyze_pathway(Normalized_Counts_Data_TMM, Sample_Data, species, pathway, collections)
})

# If you need to clean up existing results
# results_list <- lapply(results_list, function(x) {
#   x[, 1:(ncol(x) - ncol(Sample_Data))]
# })

names(results_list) <- pathways

# str(results_list[["WNT"]])
# colnames(results_list[["WNT"]])
# 
# str(Sample_Data)
# rownames(Sample_Data)
# 
# print("GSVA samples:")
# print(colnames(results_list[["WNT"]]))
# print("Sample_Data samples:")
# print(rownames(Sample_Data))
library(grid)

# Function to create a heatmap for one pathway
create_heatmap <- function(gsva_scores, pathway, sample_data, total_height = unit(8, "inch")) {
  # Ensure sample_data is in the same order as gsva_scores columns
  sample_data <- sample_data[colnames(gsva_scores), ]
  
  # Create color vectors for each annotation
  # slide_colors <- setNames(rainbow(length(unique(sample_data$SlideName))), 
  #                          unique(sample_data$SlideName))
  slide_colors <- setNames(c("#FF3158", "#42B858", "#4268F4"), unique(sample_data$SlideName))
  neuron_colors <- setNames(c("#707070", "#F656F4"), unique(sample_data$Neuron))
  
  # Create top annotation
  ha_top <- HeatmapAnnotation(
    SlideName = anno_simple(sample_data$SlideName, col = slide_colors),
    Neuron = anno_simple(sample_data$Neuron, col = neuron_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    gap = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(c(4, 4), "mm")
  )
  
  # Adjust row names if they're too long
  row_names <- rownames(gsva_scores)
  if(max(nchar(row_names)) > 50) {
    row_names <- substr(row_names, 1, 50)
  }
  
  # Calculate total width
  n_cols <- ncol(gsva_scores)
  total_width <- unit(max(15, n_cols * 0.5), "cm")  # Adjust the multiplier (0.5) as needed
  
  # Create main heatmap
  ht <- Heatmap(gsva_scores,
                name = paste(pathway, "GSVA Scores"),
                column_title = pathway,
                row_names_gp = gpar(fontsize = 8),
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                row_names_max_width = unit(10, "cm"),  # Increased width for row names
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 6),
                top_annotation = ha_top,
                column_names_side = "bottom",
                column_names_rot = 90,
                column_split = sample_data$SlideName,
                column_gap = unit(2, "mm"),
                show_column_dend = TRUE,
                show_row_dend = FALSE,
                width = total_width,  # Set total width
                heatmap_legend_param = list(
                  title = "Scores",
                  direction = "vertical",
                  title_position = "topcenter",
                  legend_height = unit(4, "cm")
                )
  )
  
  # Create annotation legends
  annotation_legend <- packLegend(
    Legend(title = "SlideName", at = names(slide_colors), legend_gp = gpar(fill = slide_colors)),
    Legend(title = "Neuron", at = names(neuron_colors), legend_gp = gpar(fill = neuron_colors)),
    direction = "vertical",
    gap = unit(1, "cm")
  )
  
  # Draw the heatmap with all legends on the left
  draw(ht, 
       heatmap_legend_side = "left",
       annotation_legend_side = "left",
       annotation_legend_list = annotation_legend,
       height = total_height,
       padding = unit(c(2, 20, 2, 10), "mm"),  # top, right, bottom, left padding
       merge_legend = FALSE)
}


# Define the SlideName variables you want to keep
selected_slides <- c("BFP", "FP", "Sympa")  # Replace with your desired slide names

# Create and draw heatmaps for each pathway
for (pathway in pathways) {
  gsva_scores <- results_list[[pathway]]
  
  # Subset Sample_Data and gsva_scores
  sample_data_subset <- Sample_Data[Sample_Data$SlideName %in% selected_slides, ]
  common_samples <- intersect(colnames(gsva_scores), rownames(sample_data_subset))
  
  if (length(common_samples) == 0) {
    warning(paste("No common samples found for pathway:", pathway))
    next
  }
  
  gsva_scores <- gsva_scores[, common_samples]
  sample_data_subset <- sample_data_subset[common_samples, ]
  
  # Calculate total height based on number of rows
  num_rows <- nrow(gsva_scores)
  total_height <- unit(min(16, max(10, num_rows * 0.3)), "inch")  # Adjusted height range
  
  pdf(paste0(pathway, "_TMM.pdf"), width = 28, height = 16)  # Increased width and height
  create_heatmap(gsva_scores, pathway, sample_data_subset, total_height)
  dev.off()
  print(paste("Heatmap for", pathway, "pathway with selected slides saved."))
}





## Volcano Plot of GeneSets ----
results_list$WNT
MouseHallmark <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

MouseHallmarkEnrichment <- GSEA(gene_list, 
                                TERM2GENE = MouseHallmark,
                                minGSSize = 10,
                                maxGSSize = 500,
                                eps = 1e-10,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH",
                                verbose = TRUE)

MouseHallmarkEnrichment_Results <- as.data.frame(MouseHallmarkEnrichment)
MouseHallmarkEnrichment_Results[,"SIG"] <- ifelse(MouseHallmarkEnrichment_Results[,"p.adjust"]<=0.05,as.numeric(sign(MouseHallmarkEnrichment_Results[,"NES"])),NA)
MouseHallmarkEnrichment_Results <- MouseHallmarkEnrichment_Results[order(MouseHallmarkEnrichment_Results[,"NES"],decreasing = T),]
MouseHallmarkEnrichment_Results <- MouseHallmarkEnrichment_Results[order(MouseHallmarkEnrichment_Results[,"SIG"],decreasing = T),]
# Convert the list column to character strings
#GO_ENRICHMENT$leadingEdge <- sapply(GO_ENRICHMENT$leadingEdge, function(x) paste(x, collapse = ","))
write.table(MouseHallmarkEnrichment_Results,paste(DOWNSTREAM_ANALYSIS_DIR,"GSEA_MouseHallmarkEnrichment.txt",sep=""),sep="\t",row.names = F,col.names = T,quote = T)




#### Delete: Normalization Definitions -----

Q3 (Upper Quartile) Normalization:
  
  Q3 normalization, also known as Upper Quartile normalization, is a method used to adjust for differences in sequencing depth between samples. Here's how it works:

    For each sample, remove genes with zero counts across all samples.
    Calculate the 75th percentile (upper quartile) of the remaining counts for each sample.
    Divide all counts in a sample by that sample's upper quartile value.
Multiply by the mean upper quartile across all samples to maintain the original scale.

This method is less sensitive to extreme values than total count normalization and can perform better when a large proportion of genes are not expressed in one of the conditions being compared.

TMM (Trimmed Mean of M-values) Normalization:
  
  TMM normalization is based on the assumption that most genes are not differentially expressed. It calculates scaling factors for each sample that can be used to adjust library sizes. Here's a simplified explanation of how it works:

    Choose a reference sample (often the sample whose library size is closest to the mean).
    For each sample, calculate log-fold-changes (M-values) and absolute expression levels (A-values) relative to the reference sample.
    Remove the genes with the highest and lowest M-values (typically 30%) and those with the highest and lowest A-values (typically 5%).
    Calculate a weighted mean of the remaining M-values, where the weights are inversely proportional to the approximate asymptotic variances.
    Convert this mean back to a scaling factor for the original library size.

TMM is designed to be robust against the presence of highly expressed genes and can perform well in situations where there are significant differences in RNA composition between samples. Both methods aim to make samples more comparable by adjusting for technical differences in sequencing depth or library composition. The choice between them (or other methods) often depends on the specific characteristics of your data and experimental design. TMM is generally considered more robust and is often the default in popular RNA-seq analysis packages like edgeR, while Upper Quartile normalization is simpler and can be effective in many cases.