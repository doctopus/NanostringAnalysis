library(ggradar)
library(scales)
library(dplyr)
library(tidyr)
library(withr)

# Prepare the data
data <- as.data.frame(seo_qc@colData)

# Define QC parameters and thresholds
QCParams <- data.frame(
  numeric_cols = c("AOISurfaceArea", "AOINucleiCount", "SequencingSaturation", "RawReads", "lib_size", "AlignedReads"),
  thresholds = c(25000, 300, 50, 2.5e+06, 2e+05, 2500000)
)

# Select relevant data (numeric columns and factors for coloring)
plot_data <- data %>%
  dplyr::select(all_of(c(QCParams$numeric_cols, "SlideName", "ROILabel", "Neuron", "Tissue")))

# Normalize the data and calculate normalized thresholds
normalize_with_threshold <- function(x, threshold) {
  all_values <- c(x, threshold)
  (all_values - min(all_values, na.rm = TRUE)) / (max(all_values, na.rm = TRUE) - min(all_values, na.rm = TRUE))
}

normalized_thresholds <- numeric(length(QCParams$numeric_cols))
for (i in seq_along(QCParams$numeric_cols)) {
  col <- QCParams$numeric_cols[i]
  threshold <- QCParams$thresholds[i]
  normalized_values <- normalize_with_threshold(plot_data[[col]], threshold)
  plot_data[[col]] <- normalized_values[1:nrow(plot_data)]
  normalized_thresholds[i] <- normalized_values[length(normalized_values)]
}

# Ensure that there are no missing or infinite values in numeric columns
plot_data[QCParams$numeric_cols] <- lapply(plot_data[QCParams$numeric_cols], function(x) {
  x[is.na(x) | !is.finite(x)] <- 0
  return(x)
})

# Convert the data to a long format for ggradar
plot_data_long <- plot_data %>%
  pivot_longer(cols = all_of(QCParams$numeric_cols), names_to = "Metric", values_to = "Value") %>%
  mutate(Type = "Data")

# Prepare the threshold data
threshold_data <- data.frame(
  Metric = QCParams$numeric_cols,
  Value = normalized_thresholds,
  Type = "Threshold"
)

# Combine the threshold and plot data
plot_data_combined <- bind_rows(plot_data_long, threshold_data)

# Convert to wide format for ggradar, only keeping numeric columns for the radar plot
plot_data_wide <- plot_data_combined %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  dplyr::select(-Type)

# Ensure no missing or infinite values in numeric columns after pivot
plot_data_wide[QCParams$numeric_cols] <- lapply(plot_data_wide[QCParams$numeric_cols], function(x) {
  x[is.na(x) | !is.finite(x)] <- 0
  return(x)
})

# Add the max and min values to the numeric columns only
plot_data_wide_numeric <- plot_data_wide %>%
  dplyr::select(all_of(QCParams$numeric_cols))

plot_data_wide_final <- rbind(
  rep(1, ncol(plot_data_wide_numeric)),  # max values for each metric
  rep(0, ncol(plot_data_wide_numeric)),  # min values for each metric
  plot_data_wide_numeric
)

# Reattach the first column (SlideName) for plotting
plot_data_wide_final <- cbind(
  "Group" = c("Max", "Min", plot_data_wide$SlideName[-c(1:2)]),  # Align SlideName with the numeric data
  plot_data_wide_numeric
)


# Generate the radar plot
ggradar(plot_data_wide_final,
        values.radar = c("0%", "50%", "100%"),
        group.colours = c("red", "forestgreen"),
        grid.line.width = 0.5,
        axis.label.size = 4,
        group.line.width = 1,
        fill = TRUE,  # Fill the area for the threshold
        fill.alpha = 0.3,  # Transparency for the threshold area
        legend.title = "Legend") +
  geom_line(aes(x = `Group`, y = value, color = `Group`), linewidth = 1)  # Use linewidth instead of size


###-----
###GGPLot2 version
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming plot_data_wide_final is already prepared and has a Group column

# Reshape the data for plotting with ggplot2
plot_data_long <- plot_data_wide_final %>%
  pivot_longer(-Group, names_to = "Metric", values_to = "Value")

# Manually create the radar plot using ggplot2
ggplot(plot_data_long, aes(x = Metric, y = Value, group = Group, color = Group)) +
  geom_polygon(aes(fill = Group), alpha = 0.3) +  # Fill the polygon area with transparency
  geom_line(size = 1) +  # Draw the lines connecting the points
  coord_polar() +  # Convert to polar coordinates for the radar chart
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(size = 12),  # Adjust text size for clarity
        legend.position = "bottom") +  # Place legend at the bottom
  labs(title = "QC Parameters with Thresholds on a Normalized Scale", 
       fill = "Group", color = "Group")  # Add titles and labels

#### test data
library(ggradar)
library(dplyr)

# Create a simplified test dataset
test_data <- data.frame(
  Group = c("Max", "Min", "Sample1", "Sample2", "Sample2"),
  Metric1 = c(1, 0, 0.8, 0.6, 0.4),
  Metric2 = c(1, 0, 0.7, 0.4, 0.5),
  Metric3 = c(1, 0, 0.9, 0.7, 0.6),
  Metric4 = c(1, 0, 0.85, 0.65, 0.7),
  Metric5 = c(1, 0, 0.75, 0.5, 0.8)
)

# Generate the radar plot
ggradar(test_data,
        values.radar = c("0%", "50%", "100%"),
        group.colours = c("red", "forestgreen", "blue", "orange"),
        grid.line.width = 0.5,
        axis.label.size = 4,
        group.line.width = 1,
        fill = TRUE,  # Fill the area for the threshold
        fill.alpha = 0.3,  # Transparency for the threshold area
        legend.title = "Legend")

###=====
library(ggradar)
library(dplyr)
library(tidyr)
 
# Ensure Group column is unique for each row by using ROILabel
plot_data_wide_final <- plot_data_wide_final %>%
  mutate(Group = ROILabel)  # Use ROILabel for unique identification

# Subset data to include the threshold row separately
threshold_row <- plot_data_wide_final %>%
  filter(Group == "Threshold")

# Remove the threshold from the main data
plot_data_main <- plot_data_wide_final %>%
  filter(Group != "Threshold")

# Generate the radar plot for the threshold area first
ggradar(threshold_row,
        values.radar = c("0%", "50%", "100%"),
        group.colours = "red",  # Use a single color for threshold
        grid.line.width = 0.5,
        axis.label.size = 4,
        group.line.width = 2,
        fill = TRUE,  # Fill the area for the threshold
        fill.alpha = 0.3,  # Transparency for the threshold area
        legend.title = "Legend")

# Overlay the lines for the main data
ggradar(plot_data_main,
        values.radar = c("0%", "50%", "100%"),
        group.colours = scale_colour_manual(values = c("orange", "blue", "green")),  # Customize colors
        grid.line.width = 0.5,
        axis.label.size = 4,
        group.line.width = 1,
        fill = FALSE,  # No fill for the other areas
        add = TRUE,  # Overlay on the existing plot
        legend.title = "Legend")

##_____
Goal Recap:
  
  Unique Observations: Each row represents a unique observation, and you want to plot each one with a different color.
Area Highlight: Only the area formed by the threshold values should be painted with transparency, while other observations are connected by lines.
Handling Non-Unique Group Values: Since Group has duplicates, we need to ensure each row is treated uniquely, but colored according to another factor (e.g., Neuron or Tissue).

Solution Approach:
  
  Assign Unique Identifiers: Use the ROILabel column as a unique identifier for each row.
Plot with ggradar:
  Use ROILabel for unique observation identification.
Color the lines based on one of your factor columns (e.g., Neuron or Tissue).
Highlight only the threshold area with transparency.

Code Implementation:
  
  Here’s how you can achieve this with ggradar:
  
  Explanation:
  
  Unique Grouping:
  The Group column is redefined to use ROILabel as the unique identifier for each row.

Separate Threshold Plotting:
  The threshold_row data is plotted first with ggradar and filled with transparency to highlight the threshold area.

Overlaying Lines:
  The main data (plot_data_main) is plotted with fill = FALSE to avoid filling the area, and the add = TRUE parameter is used to overlay this on top of the threshold plot.

Custom Colors:
  Customize the line colors using scale_colour_manual() to ensure each observation is colored according to your factor of interest (e.g., Neuron).

Additional Notes:
  
  Plot Overlays: The add = TRUE parameter is used to overlay the plot, but if this doesn’t work as expected with ggradar, consider plotting the threshold and main data separately in two steps, with one radar plot drawn after the other using the plot() function to maintain layering.
Color Customization: Adjust the colors in scale_colour_manual() based on the number of factors and desired color palette.

This approach should meet your requirements for correctly plotting each observation with different colors, highlighting the threshold area with transparency, and ensuring each observation is treated uniquely.