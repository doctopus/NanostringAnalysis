library(fmsb)
library(dplyr)

# Load data from the correct source
data <- as.data.frame(seo_qc@colData)

# Define QC parameters and thresholds
QCParams <- data.frame(
  numeric_cols = c("AOISurfaceArea", "AOINucleiCount", "SequencingSaturation", "RawReads", "lib_size", "AlignedReads"),
  thresholds = c(25000, 300, 50, 2.5e+06, 2e+05, 2500000)
)

# Prepare the data for the spider plot
plot_data <- data %>%
  dplyr::select(all_of(c(QCParams$numeric_cols, "SlideName", "ROILabel", "Neuron", "Tissue")))

# Function to normalize data including threshold
normalize_with_threshold <- function(x, threshold) {
  all_values <- c(x, threshold)
  (all_values - min(all_values, na.rm = TRUE)) / (max(all_values, na.rm = TRUE) - min(all_values, na.rm = TRUE))
}

# Normalize the data and calculate normalized thresholds
normalized_thresholds <- numeric(length(QCParams$numeric_cols))
for (i in seq_along(QCParams$numeric_cols)) {
  col <- QCParams$numeric_cols[i]
  threshold <- QCParams$thresholds[i]
  normalized_values <- normalize_with_threshold(plot_data[[col]], threshold)
  plot_data[[col]] <- normalized_values[1:nrow(plot_data)]
  normalized_thresholds[i] <- normalized_values[length(normalized_values)]
}

# Create max and min rows with correct column names
max_values <- setNames(rep(1, length(QCParams$numeric_cols)), QCParams$numeric_cols)
min_values <- setNames(rep(0, length(QCParams$numeric_cols)), QCParams$numeric_cols)
max_min <- rbind(max_values, min_values)

# Add threshold row
threshold_row <- setNames(normalized_thresholds, QCParams$numeric_cols)

# Combine max_min with threshold and plot_data
plot_data_final <- rbind(max_min, threshold_row, plot_data[QCParams$numeric_cols])

# Remove any rows with NA or NaN values
plot_data_final <- na.omit(plot_data_final)

# Ensure all values are numeric
plot_data_final[] <- lapply(plot_data_final, as.numeric)

# Set colors manually
color_map <- c("NF-H-" = "darkorange", "NF-H+" = "purple")

# Create the spider plot
par(mar = c(2, 2, 2, 2))
radarchart(
  plot_data_final,
  pcol = c("forestgreen", color_map[plot_data$Neuron]),
  pfcol = c(scales::alpha("forestgreen", 0.3), rep(NA, nrow(plot_data))),
  plwd = c(2, rep(1, nrow(plot_data))),
  plty = c(1, rep(1, nrow(plot_data))),
  cglcol = "grey",
  cglty = 1,
  axistype = 1,
  caxislabels = seq(0, 1, 0.25),
  cglwd = 1,
  vlcex = 0.8,
  axislabcol = "black",
  seg = 5
)

radarchart(
  plot_data_final,
  cglcol = "grey"
  # pcol = c("red", color_map[plot_data$Neuron]),
  # pfcol = c(scales::alpha("red", 0.3), rep(NA, nrow(plot_data)))
)
# Add a legend
legend(
  # "bottomleft",
  # 2.8, 19,
  x=c(0.8, 1.04),
  y=c(0.2, 0.3),
  # text.width = 15,
  legend = c("Threshold", names(color_map)),
  col = c("red", color_map),
  lwd = c(2, 1, 1),
  bty = "n"
)

# Add title
title("QC Parameters with Thresholds on a Normalized Scale")