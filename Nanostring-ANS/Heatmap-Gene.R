###Plot HandPicked Gene Heatmap from Normalized Counts
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr) #For %>% Operator

create_heatmap <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), ]
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#FF3158", "#42B858", "#4268F4"), unique(sample_data$SlideName)) #, "#C9F52A", "#A0582A"
  neuron_colors <- setNames(c("#707070", "#F656F4"), unique(sample_data$Neuron))
  
  # FOR z-Score Normalized data: Create a diverging color palette
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  # Set breaks for the color scale
  max_abs <- max(abs(count_scores))
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  
  # Create top annotation with custom colors
  ha_top <- HeatmapAnnotation(
    SlideName = anno_simple(sample_data$SlideName, col = slide_colors),
    Neuron = anno_simple(sample_data$Neuron, col = neuron_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )
  # Create main heatmap
  ht <- Heatmap(count_scores,
                name = paste(pathway, "Scores"),
                col = colorRamp2(breaks, colors), #Added for scaled_data; not needed otherwise
                column_title = pathway,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 12),  # Adjust font size of gene names
                top_annotation = ha_top,
                column_split = sample_data$SlideName,
                heatmap_legend_param = list(title = "Expression")
  )
  # Create annotation legends
  annotation_legend <- packLegend(
    Legend(title = "SlideName", at = names(slide_colors), legend_gp = gpar(fill = slide_colors)),
    Legend(title = "Neuron", at = names(neuron_colors), legend_gp = gpar(fill = neuron_colors)),
    direction = "vertical",
    gap = unit(1, "cm")
  )
  num_rows <- nrow(count_scores)
  total_height <- unit(min(15, max(10, num_rows * 0.4)), "inch")
  # Draw the heatmap with all legends on the left
  
  draw(ht, 
       heatmap_legend_side = "left",
       annotation_legend_side = "left",
       annotation_legend_list = annotation_legend,
       height = total_height,
       padding = unit(c(2, 20, 2, 10), "mm"),  # top, right, bottom, left padding
       merge_legend = FALSE)
}

plot_and_save_heatmap <- function(normalizedCountsData, sample_data_subset, pathway, output_file) {
  # Calculate height based on number of genes
  # height <- max(12, nrow(normalizedCountsData) * 0.2)  # Adjust the multiplier (0.2) as needed
  height <- 18
  width <- 16
  
  pdf(output_file, width = width, height = height)  # Increased width to accommodate legends
  create_heatmap(normalizedCountsData, pathway, sample_data_subset)
  dev.off()
  print(paste("Heatmap for", pathway, "pathway saved to", output_file))
}

#Load Input Data
custom_genes <- c("Fzd1", "Fzd2", "Fzd3", "Fzd4", "Fzd5", "Fzd6", "Fzd7", "Fzd8", "Fzd9", "Fzd10","Lrp5", "Lrp6", "Rnf43", "Znrf3", "Lgr5", 
                  "Celsr1", "Celsr2", "Celsr3", "Vangl1", "Vangl2", "Ror1", "Ror2", "Ryk", "Ptk7", "Musk", "Sdc1", "Gpc1", "Orai1", "Orai3", "Ctbp1", "Ctbp2", "Pygo1", "Pygo2",  
                  "Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt4", "Wnt5b", "Wnt6", "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b", "Wnt11", "Wnt16",
                  "Tmed2", "Gpr177", "Notum", "Norrin", "Cer1", "Sfrp1",  "Sfrp2", "Sfrp5", "Wif", "Sost", "Dkk1", "Dkk4", "Draxin", "Igfbp4", "Rspo3", 
                  "Snai2", "Sox2", "Sox17", "Adamts5", "Adam11", "Ccnd1", "Ccna2", "Vegfa", "Npy")

custom_genes <- c("Nlgn3", "Nlgn1", "Npy", "Nlgn2")

custom_genes <- c("Mmp2", "Jun", "Adam17", "Adam10", "Nlgn2", "Mmp9", "Nlgn3", "Adrb2", "Adrb1", "Npy1r", "Nlgn1", "Npy", "Nrxn1", "Nrxn2", "Nrxn3", "Npy5r", "Slc18a3", "Npy2r", "Fos")

custom_genes <- c("Fzd10","Wnt3","Wnt6","Wnt5b","Wnt9b","Wnt11","Rspo3","Dkk4", "Draxin", "Ngf","Snai2","Sox2","Sox17","Adamts5","Adam11")

custom_genes <- c("Zeb1", "Tcf7", "Snail", "Myc", "Abc", "Ccnd1", "Cdh2", "Fak", "Akt", "Jnk", "Erk", "Mapk", "Stat3", "Axin2", "Rspo3",
                  "Wnt11", "Wnt5b", "Wnt3", "Fzd10")
#Diagnose duplicates
# gene_names <- data.frame(rownames=rownames(Normalized_Counts_Data_TMM))


normalizedCountsData <- Normalized_Counts_Data_TMM #[INPUT_NEEDED] Or Normalized_Counts_Data_TMM

normalizedCountsData <- normalizedCountsData[custom_genes,]
normalizedCountsData <- normalizedCountsData[!grepl("^NA", rownames(normalizedCountsData)),]


#Per Gene Scaling of filtered Normalized count ----
# Function to scale each row (gene) (Decide if to scale across all samples or after the filtering below)
scale_rows <- function(x) {
  (x - mean(x)) / sd(x)
}
# Apply scaling to each row
scaled_data <- t(apply(normalizedCountsData, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data[is.nan(scaled_data)] <- 0

#Renamed back again to normalizedCountsData for simplicity of final code
normalizedCountsData <- scaled_data


#Process Sample Data----
sampleData <- Sample_Data #[INPUT_NEEDED]
sampleData <- sampleData %>% mutate(SlideName = gsub("BBP", "BFP", SlideName))
# Define the SlideName variables you want to keep
selected_slides <- c("BFP", "FP", "Sympa")#, "IBTP", "SP")  # Replace with your desired slide names

# Subset sampleData and normalizedCountsData
sample_data_subset <- sampleData[sampleData$SlideName %in% selected_slides, ]
common_samples <- intersect(colnames(normalizedCountsData), rownames(sample_data_subset))

if (length(common_samples) == 0) {
  stop("No common samples found")
}

normalizedCountsData <- normalizedCountsData[, common_samples]

normalizedCountsData <- as.matrix(normalizedCountsData)
sample_data_subset <- sample_data_subset[common_samples, ]

# Plot and save the heatmap


plot_and_save_heatmap(normalizedCountsData, 
                      sample_data_subset, 
                      "Wnt Pathway Gene Expression", 
                      "Custom_Pathway_TMM.pdf")
