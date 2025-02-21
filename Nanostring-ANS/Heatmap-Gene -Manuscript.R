###Plot HandPicked Gene Heatmap from Normalized Counts -All legends in one column
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr) #For %>% Operator

process_expression_data <- function(raw_counts, custom_genes, sample_data, selected_slides) {
  # Filter genes
  filtered_counts <- raw_counts[custom_genes,]
  filtered_counts <- filtered_counts[!grepl("^NA", rownames(filtered_counts)),]
  
  # Scale data
  scaled_counts <- scale_expression_data(filtered_counts)
  
  # Process sample data
  processed_sample_data <- process_sample_metadata(sample_data, selected_slides)
  
  # Ensure data alignment
  common_samples <- intersect(colnames(scaled_counts), rownames(processed_sample_data))
  if (length(common_samples) == 0) {
    stop("No common samples found between count data and sample metadata")
  }
  
  aligned_counts <- as.matrix(scaled_counts[, common_samples])
  aligned_sample_data <- processed_sample_data[common_samples, ]
  
  return(list(
    counts = aligned_counts,
    metadata = aligned_sample_data
  ))
}

scale_expression_data <- function(count_data) {
  scale_rows <- function(x) {
    (x - mean(x)) / sd(x)
  }
  scaled_data <- t(apply(count_data, 1, scale_rows))
  scaled_data[is.nan(scaled_data)] <- 0
  return(scaled_data)
}

process_sample_metadata <- function(sample_data, selected_slides) {
  sample_data %>% 
    mutate(SlideName = gsub("BBP", "BFP", SlideName)) %>%
    filter(SlideName %in% selected_slides)
}

create_heatmap <- function(count_scores, pathway, sample_data) {
  # Prepare annotations and colors
  annotations <- prepare_annotations(count_scores, sample_data)
  color_scheme <- create_color_scheme(count_scores, sample_data)
  
  # Create heatmap object
  heatmap_obj <- create_main_heatmap(count_scores, pathway, sample_data, 
                                     annotations$ha_top, color_scheme)
  
  # Create and arrange legends with all necessary parameters
  legends <- create_legends(
    color_scheme$slide_colors, 
    color_scheme$neuron_colors,
    count_scores,
    color_scheme
  )
  
  # Draw final heatmap
  draw_complete_heatmap(heatmap_obj, legends, count_scores)
}

prepare_annotations <- function(count_scores, sample_data) {
  sample_data <- sample_data[colnames(count_scores), ]
  
  ha_top <- HeatmapAnnotation(
    SlideName = anno_simple(sample_data$SlideName, 
                            col = setNames(c("#FF3158", "#42B858", "#4268F4"), 
                                           unique(sample_data$SlideName))),
    Neuron = anno_simple(sample_data$Neuron, 
                         col = setNames(c("#F656F4", "#707070"), 
                                        unique(sample_data$Neuron))),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 16, fontface = "bold") #  Adjust font size and style
  )
  
  return(list(ha_top = ha_top))
}

create_color_scheme <- function(count_scores, sample_data) {
  slide_colors <- setNames(c("#FF3158", "#42B858", "#4268F4"), 
                           unique(sample_data$SlideName))
  neuron_colors <- setNames(c("#F656F4", "#707070"), 
                            unique(sample_data$Neuron))
  
  expression_colors <- colorRampPalette(c("blue", "white", "red"))(101)
  max_abs <- max(abs(count_scores))
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  return(list(
    slide_colors = slide_colors,
    neuron_colors = neuron_colors,
    expression_colors = expression_colors,
    breaks = breaks
  ))
}

create_main_heatmap <- function(count_scores, pathway, sample_data, ha_top, color_scheme) {
  Heatmap(count_scores,
          name = paste(pathway, "Scores"),
          col = colorRamp2(color_scheme$breaks, color_scheme$expression_colors),
          # column_title = pathway,
          column_title = c("Before First Peak", "First Peak", "Sympathectomy"),
          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
          
          # Split column settings
          column_split = sample_data$SlideName,
          column_gap = unit(2, "mm"),
          column_title_rot = 0,
          column_names_gp = gpar(fontsize = 18, fontface = "bold"),  # This will style the split labels
          
          # Other settings
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = FALSE,  # Hide individual column names
          row_names_gp = gpar(fontsize = 16, fontface = "bold"),
          top_annotation = ha_top,
          show_heatmap_legend = FALSE
  )
}

create_legends <- function(slide_colors, neuron_colors, count_scores, color_scheme) {
  # Create expression legend
  expression_legend <- Legend(
    title = "Expression",
    col_fun = colorRamp2(color_scheme$breaks, color_scheme$expression_colors),
    at = c(-max(abs(count_scores)), 0, max(abs(count_scores))),
    labels = c("Low", "Medium", "High"),
    title_gp = gpar(fontsize = 20, fontface = "bold"),  # Increase title size
    labels_gp = gpar(fontsize = 18),  # Increase labels size
    title_gap = unit(8, "mm"),  # Add gap between title and scale
    row_gap = unit(5, "mm")     # Add gap between items in the legend
  )
  # Define the desired order for slides
  slide_order <- c("BFP", "FP", "Sympa")  # Replace with your desired order
  
  # Create slide legend with specified order
  slide_legend <- Legend(
    title = "SlideName", 
    at = slide_order,  # Use the ordered names
    legend_gp = gpar(fill = slide_colors[slide_order]),  # Reorder colors to match
    title_gp = gpar(fontsize = 20, fontface = "bold"),
    labels_gp = gpar(fontsize = 18),
    title_gap = unit(8, "mm"),
    row_gap = unit(5, "mm")
  )
  
  neuron_legend <- Legend(
    title = "Neuron", 
    at = names(neuron_colors), 
    legend_gp = gpar(fill = neuron_colors),
    title_gp = gpar(fontsize = 20, fontface = "bold"),  # Increase title size
    labels_gp = gpar(fontsize = 18),  # Increase labels size
    title_gap = unit(8, "mm"),  # Add gap between title and items
    row_gap = unit(5, "mm")     # Add gap between items
  )
  
  # Pack all legends vertically
  packLegend(
    expression_legend,
    slide_legend,
    neuron_legend,
    direction = "vertical",
    gap = unit(2, "cm")
  )
}

draw_complete_heatmap <- function(heatmap_obj, legends, count_scores) {
  total_height <- unit(min(18, max(10, nrow(count_scores) * 0.3)), "inch")
  
  draw(heatmap_obj,
       heatmap_legend_side = "left",
       annotation_legend_side = "left",  # Changed from separate legend sides
       annotation_legend_list = legends,
       height = total_height,
       padding = unit(c(2, 20, 2, 10), "mm"),
       merge_legend = TRUE)  # Changed to TRUE to ensure legends are merged
}


save_heatmap <- function(count_data, sample_data, pathway, output_file, 
                         width = 16, height = 18) {
  # Create heatmap object
  heatmap_obj <- create_heatmap(count_data, pathway, sample_data)
  
  # Save to PDF
  pdf(output_file, width = width, height = height)
  draw(heatmap_obj)
  dev.off()
  
  # Display in RStudio viewport
  print(heatmap_obj)
  
  message(sprintf("Heatmap for %s pathway saved to %s", pathway, output_file))
}

# Usage:
processed_data <- process_expression_data(
  Normalized_Counts_Data_TMM,
  custom_genes,
  Sample_Data,
  c("BFP", "FP", "Sympa")
)
#
save_heatmap(
  processed_data$counts,
  processed_data$metadata,
  "Wnt Pathway Gene Expression",
  "20250212-29 Custom_Pathway_TMM.pdf"
)
