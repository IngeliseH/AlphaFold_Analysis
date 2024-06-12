library(tidyverse)

# Function to load and process data
process_metrics_data <- function(file_path) {
  data <- read_csv(file_path)
  data <- data %>%
    mutate(
      ipTM = as.numeric(ipTM),
      min_PAE = as.numeric(min_PAE),
      pdockq = as.numeric(pDockQ)
    )
  return(data)
}

# Function to create a single scatter plot with flexible inputs and optional coloring
plot_data <- function(data, x_var, y_var, color_var = NULL, file_name = "plot.png") {
  # Check and create the directory if it does not exist
  plot_dir <- dirname(file_name)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Define the path and settings for the PNG output
  png(file_name, width = 2000, height = 2000, res = 300)  # Set to 800x800 pixels, 300 DPI

  # Start the ggplot
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.2, aes_string(color = color_var)) +  # Conditional aes_string for color
    labs(title = paste("Relationship between", x_var, "and", y_var),
         x = x_var, y = y_var) +
    theme_minimal() +
    theme(text = element_text(size = 12),  # Increase global text size
      axis.title = element_text(size = 12),  # Increase axis title size
      axis.text = element_text(size = 11),  # Increase axis text size
      plot.title = element_text(size = 14, face = "bold")
    )  # Bold and larger title

  # Adjusting plot dimensions to be square
  p <- p + theme(aspect.ratio = 1)

  # Check if color_var is not NULL and a valid column name
  if (!is.null(color_var) && color_var %in% names(data)) {
    if (is.factor(data[[color_var]])) {
      p <- p + scale_color_viridis_d()
    } else {
      p <- p + scale_color_viridis_c()
    }
  } else {
    p <- p + scale_color_manual(values = "black")  # Default to blue if no valid color_var
  }

  # Print the plot
  print(p)

  # Close the PNG device to finalize and save the plot
  dev.off()
}

file_path <- "alphafold_predictions_results.csv"
data <- process_metrics_data(file_path)
summary(data)
# Plot min_PAE vs iptm
plot_data(data, "min_PAE", "ipTM", file_name = "plots/min_PAE_vs_iptm.png")
# Plot min_PAE vs pdockq
plot_data(data, "min_PAE", "pdockq", file_name = "plots/min_PAE_vs_pdockq.png")
# Plot pdockq vs iptm
plot_data(data, "pdockq", "ipTM", file_name = "plots/pdockq_vs_iptm.png")
# Plot pdockq vs iptm with pdockq as color
plot_data(data, "min_PAE", "ipTM", "pdockq", file_name = "plots/min_PAE_vs_iptm_vs_pdockq_color.png")
