#NOT WORKING
#Failure in merging or in processing evidence - too much data marked as having evidence, and evidence type not working
# read lintr file
read.dcf(".lintr")

# Load necessary libraries
library(dplyr)
# For 3D plots
library(rgl)

# Function to load data supporting CSV and Excel formats with error handling
load_data <- function(file_path) {
  if (grepl("\\.csv$", file_path)) {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else if (grepl("\\.xlsx$", file_path)) {
    return(readxl::read_excel(file_path))
  } else {
    stop("Unsupported file format: ", file_path)
  }
}

# Function to convert specified columns to numeric and replace NAs with zeros
convert_to_numeric <- function(data, columns) {
  for (column in columns) {
    data[[column]] <- as.numeric(data[[column]])
  }
  return(data)
}

# Function to discretize a continuous variable into specified bins as an ordered factor
discretize_variable <- function(data, column, breaks, labels, right = TRUE) {
  data[[paste0(column, "_discrete")]] <- cut(data[[column]], breaks = breaks, labels = labels, right = right, include.lowest = TRUE)
  data[[paste0(column, "_discrete")]] <- factor(data[[paste0(column, "_discrete")]], levels = labels, ordered = TRUE)
  return(data)
}

# Process alphafold data from csv for analysis
prep_alphafold_data <- function(alphafold_path) {
  # Read in alphafold data
  alphafold_data <- load_data(alphafold_path)
  # Convert selected columns to numeric
  alphafold_data <- convert_to_numeric(alphafold_data, c("ipTM", "min_PAE", "pDockQ"))
  # Add discretized columns for selected variables
  alphafold_data <- discretize_variable(alphafold_data, "pDockQ", c(0, 0.23, 0.49, 1.0), c("0-0.23", "0.23-0.49", "0.49-1.0"))
  alphafold_data <- discretize_variable(alphafold_data, "ipTM", c(0, 0.55, 0.7, 1.0), c("0-0.55", "0.55-0.7", "0.7-1.0"))
  alphafold_data <- discretize_variable(alphafold_data, "min_PAE", c(0, 2.5, 5, 10, 30), c("0-2.5", "2.5-5", "5-10", "10-30"))
}

# Process and merge AlphaFold and experimental data into a single dataset
prep_for_experimental_comparison <- function(alphafold_path, experimental_path) {
  # Read in alphafold data
  alphafold_data <- load_data(alphafold_path)
  # Make key column for merging - Protein_Pair
  alphafold_data$Protein_Pair <- apply(alphafold_data[,c("Protein1", "Protein2")], 1, function(x) paste(sort(x), collapse = "_"))
  # Convert selected columns to numeric
  alphafold_data <- convert_to_numeric(alphafold_data, c("ipTM", "min_PAE", "pDockQ"))
  # Condense to one row per protein pair
  alphafold_data = alphafold_data %>% group_by(Protein_Pair) %>% summarise(ipTM = max(ipTM), min_PAE = min(min_PAE), pDockQ = max(pDockQ), ppv = max(ppv), Num_Consistent = max(Num_Consistent))
  # Add discretized columns for selected variables
  alphafold_data <- discretize_variable(alphafold_data, "pDockQ", c(0, 0.23, 0.49, 1.0), c("0-0.23", "0.23-0.49", "0.49-1.0"))
  alphafold_data <- discretize_variable(alphafold_data, "ipTM", c(0, 0.55, 0.7, 1.0), c("0-0.55", "0.55-0.7", "0.7-1.0"))
  alphafold_data <- discretize_variable(alphafold_data, "min_PAE", c(0, 2.5, 5, 10, 30), c("0-2.5", "2.5-5", "5-10", "10-30"))
  # Copy num_consistent as ordered factor
  alphafold_data$Num_Consistent_discrete <- factor(alphafold_data$Num_Consistent, levels = c(0, 1, 2, 3, 4), ordered = TRUE)

  # Read in experimental data
  experimental_data <- load_data(experimental_path)
  # Remove duplicate rows
  experimental_data <- experimental_data[!duplicated(experimental_data), ]
  # Create Protein_Pair key column for merging, in alphabetical order
  experimental_data$Protein_Pair <- apply(experimental_data[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")], 1, function(x) paste(sort(x), collapse = "_"))

  # Merge datasets
  merged_data = merge(alphafold_data, experimental_data, by = "Protein_Pair", all = TRUE)

  # Add a binary variable for evidence of interaction
  merged_data$Evidence <- ifelse(!is.na(merged_data$`Official.Symbol.Interactor.A`), "Yes", "No")
  merged_data$Evidence <- factor(merged_data$Evidence, levels = c("Yes", "No"))
  # Add a factored variable for evidence type
  merged_data <- merged_data %>%
    mutate(Evidence_Type = case_when(
      Evidence == "No" ~ "No",
      grepl("Structure", Experimental.System) ~ "Structure",
      TRUE ~ "Other"
    ))
  merged_data$Evidence_Type <- factor(merged_data$Evidence_Type, levels = c("Structure", "Other", "No"))

  return(merged_data)
}

# Function to perform logistic regression analysis on the data
perform_logistic_regression <- function(data, feature_columns, target_column) {
  # Initial setup and data preparation
  formula <- as.formula(paste(target_column, "~", paste(feature_columns, collapse = "+")))
  model <- glm(formula, data = data, family = "binomial")
  summary(model)
}

# Visualization functions for plotting data
# Function to create and save a plot with flexible input and customization options
plot_data <- function(data, x_var, y_var, color_var = NULL, file_name, alpha = 0.5, add_regression = FALSE) {
  plot_dir <- dirname(file_name)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = alpha, aes_string(color = color_var)) +
    labs(title = paste("Relationship between", x_var, "and", y_var), x = x_var, y = y_var) +
    theme_minimal() +
    theme(text = element_text(size = 16), axis.title = element_text(size = 16), title = element_text(size = 18))


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

  if (add_regression && !is.null(color_var)) {
    p <- p + geom_smooth(method = "lm", aes(group = color_var), se = FALSE)
  } else if (add_regression) {
    p <- p + geom_smooth(method = "lm", color = "red", se = FALSE)
  }

  ggsave(file_name, plot = p, width = 10, height = 8)
}

# Function to plot box plots of metrics versus evidence
plot_box_plots <- function(data, metrics, x_column, color_var = NULL) {
  for (metric in metrics) {
    p <- ggplot(data, aes_string(x = x_column, y = metric)) +
      geom_boxplot()

    # Adding jitter with optional coloring
    if (!is.null(color_var) && color_var %in% names(data)) {
      if (is.factor(data[[color_var]])) {
        p <- p + geom_jitter(width = 0.05, height = 0, alpha = 0.5, aes_string(color = color_var)) +
                scale_color_viridis_d()
      } else {
        p <- p + geom_jitter(width = 0.05, height = 0, alpha = 0.5, aes_string(color = color_var)) +
                scale_color_viridis_c()
      }
    } else {
      p <- p + geom_jitter(width = 0.05, height = 0, alpha = 0.5, color = "black")  # Default to black if no valid color_var
    }

    p <- p + labs(title = paste("Box Plot of", metric, "by", x_column), x = x_column, y = metric) +
      theme_minimal() +
      theme(text = element_text(size = 16), axis.title = element_text(size = 16), title = element_text(size = 18))

    file_name <- paste("plots/", metric, "_vs_", x_column, ".png", sep = "")
    ggsave(file_name, plot = p, width = 12, height = 8)
  }
}

plot_3d <- function(data, x_var, y_var, z_var, color_var) {
  # Ensure the color_var is a factor for consistent coloring
  if (!is.factor(data[[color_var]])) {
    data[[color_var]] <- as.factor(data[[color_var]])
  }

  # Colors - using a palette that scales with the number of levels in the factor
  colors <- rainbow(length(levels(data[[color_var]])))

  # Plotting the 3D scatter plot
  plot3d(
    x = data[[x_var]], 
    y = data[[y_var]], 
    z = data[[z_var]], 
    col = colors[data[[color_var]]],  # Mapping colors based on the factor levels
    xlab = x_var, 
    ylab = y_var, 
    zlab = z_var,
    main = paste("3D Plot of", x_var, y_var, z_var, "colored by", color_var)
  )

  # Add a legend to the plot
  legend3d("topright", legend = levels(data[[color_var]]), col = colors, pch = 16)
}

# Function to plot discretized variables against evidence type
plot_discrete_evidence_comparison <- function(data, discrete_vars, evidence_type_column) {
  for (var in discrete_vars) {
    p <- ggplot(data, aes_string(x = var, fill = evidence_type_column)) +
      geom_bar(position = "dodge") +
      labs(title = paste("Distribution of", var, "by", evidence_type_column), x = var, y = "Count") +
      theme_minimal()
    file_name <- paste("plots/", var, "_vs_", evidence_type_column, ".png", sep = "")
    ggsave(file_name, plot = p, width = 12, height = 8)
  }
}

plot_data_density <- function(data, x_var, y_var, file_name, add_regression = FALSE) {
  plot_dir <- dirname(file_name)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Prepare the plot using geom_count to adjust point sizes based on overlap density
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_count(aes(size = ..n.., color = ..n..), show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(0.1, 10)) +
    labs(title = paste("Density of data points between", x_var, "and", y_var),
         x = x_var, y = y_var) +
    theme_minimal() +
    theme(text = element_text(size = 16), axis.title = element_text(size = 16), title = element_text(size = 18)) +
    guides(color = guide_legend(), size = guide_legend())

  # Adding regression line if requested
  if (add_regression) {
    p <- p + geom_smooth(method = "lm", color = "red", se = FALSE)
  }

  # Save the plot to the specified file
  ggsave(file_name, plot = p, width = 10, height = 8)
}

############################################################################################################
# Usage
# Load and process data
data <- prep_alphafold_data("data/alphafold_predictions_results.csv")
summary(data)
head(data)

# Perform statistical analysis
models <- perform_logistic_regression(data, c("ipTM", "min_PAE", "pDockQ"), "Evidence")
print(models)
# Assess each metric compared to each other
# filter data to remove rows where iptm == NA
data <- data[!is.na(data$ipTM), ]
cor(data[c("ipTM", "min_PAE", "pDockQ")])
# filter to remove rows where min_PAE > 15
data_minPAE_15 <- data[data$min_PAE <= 15, ]
cor(data_minPAE_15[c("ipTM", "min_PAE", "pDockQ", "Num_Consistent")])

# Generate plots
# Plot min_PAE vs iptm
plot_data(data, "min_PAE", "ipTM", file_name = "plots/min_PAE_vs_iptm.png")
# Plot min_PAE vs pdockq
plot_data(data, "min_PAE", "pdockq", file_name = "plots/min_PAE_vs_pdockq.png")
# Plot pdockq vs iptm
plot_data(data, "pdockq", "ipTM", file_name = "plots/pdockq_vs_iptm.png")
# Plot min PAE vs iptm with pdockq as color
plot_data(data, "min_PAE", "ipTM", "pdockq", file_name = "plots/min_PAE_vs_iptm_vs_pdockq_color.png")
# Plot min PAE vs iptm with Level_Consistent as color
plot_data(data, "min_PAE", "ipTM", "Level_Consistent", file_name = "plots/min_PAE_vs_iptm_vs_Level_Consistent_color.png")
# Plot min PAE vs iptm with Num_Consistent as color
plot_data(data, "min_PAE", "ipTM", "Num_Consistent", file_name = "plots/min_PAE_vs_iptm_vs_Num_Consistent_color.png")
# Plot num_consistent vs ipTM with level_consistent as color
plot_data(data, "Num_Consistent", "ipTM", "Level_Consistent", file_name = "plots/iptm_vs_num_consistent_vs_level_consistent_color.png")
# Plot num_consistent vs min_PAE with level_consistent as color
plot_data(data, "Num_Consistent", "min_PAE", "Level_Consistent", file_name = "plots/min_PAE_vs_num_consistent_vs_level_consistent_color.png")
# Plot num_consistent vs ipTM with min_PAE as color
plot_data(data, "Num_Consistent", "ipTM", "min_PAE", file_name = "plots/iptm_vs_num_consistent_vs_min_PAE_color.png")

# 3D plot of min_PAE, ipTM, and pdockq coloured by Num_Consistent
plot_3d(data, "ipTM", "min_PAE", "pdockq", "Num_Consistent")

############################################################################################################
# Comparing alphafold results to experimental data
merged_data <- prep_for_experimental_comparison("data/alphafold_predictions_results.csv", "data/Merged_remove_duplicate_all.csv")
summary(merged_data)
head(merged_data)

# Example usage of new plotting functions
metrics <- c("ipTM", "min_PAE", "pDockQ", "Num_Consistent")
discrete_vars <- c("pDockQ_discrete", "ipTM_discrete", "min_PAE_discrete", "Num_Consistent_discrete")

# Plotting box plots
plot_box_plots(merged_data, metrics, "Evidence", "ipTM")
plot_box_plots(merged_data, metrics, "Evidence_Type", "ipTM")

# Plotting discrete variable comparisons
plot_discrete_evidence_comparison(merged_data, discrete_vars, "Evidence_Type")

# save merged data to csv 'data/processed_data.csv'
write.csv(merged_data, "data/processed_data.csv", row.names = FALSE)

############################################################################################################
# comparing pae png assessed ROP to structure pae assessed ROP
data_2 <- prep_alphafold_data("data/alphafold_predictions_png_ROP.csv")
summary(data_2)
plot_data(data_2, "Num_Consistent", "png_ROP", "ipTM", "plots/Num_Consistent_vs_Num_Consistent_png_vs_ipTM_color.png", alpha = 0.01)
plot_data_density(data_2, "Num_Consistent", "png_ROP", "plots/Num_Consistent_vs_Num_Consistent_png_density.png", add_regression = TRUE)
