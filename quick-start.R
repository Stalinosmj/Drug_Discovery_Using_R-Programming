# Quick Start Script for Drug Discovery ML Project
# Run this script to set up and execute the entire pipeline

# ============================================================================
# SETUP AND CONFIGURATION
# ============================================================================

cat("
╔════════════════════════════════════════════════════════════════╗
║   Bioinformatics Drug Discovery - Quick Start                 ║
║   Machine Learning Pipeline for SARS-CoV-2 Inhibitors        ║
╚════════════════════════════════════════════════════════════════╝
\n")

# Clear environment
rm(list = ls())
gc()

# Check R version
if(getRversion() < "4.3.0") {
  stop("R version 4.3.0 or higher required. Current version: ", getRversion())
}
#Install Java 21
# In ~/.Rprofile or Rprofile.site
Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jdk-21")
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Program Files/Java/jdk-21/bin/server", sep = ";"))


# ============================================================================
# INSTALL REQUIRED PACKAGES (First Time Only)
# ============================================================================

cat("\n📦 Checking and installing required packages...\n")

# Function to install if needed
install_if_needed <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages) > 0) {
    cat("Installing:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, dependencies = TRUE, repos = "https://cloud.r-project.org/")
  } else {
    cat("✓ All CRAN packages already installed\n")
  }
}

# Core packages
core_packages <- c(
  "tidyverse", "data.table", "janitor",
  "rcdk", "rcdklibs", "webchem",
  "ranger", "xgboost", "caret",
  "future", "future.apply", "furrr", "progressr", "doParallel",
  "memuse", "pryr",
  "httr", "jsonlite", "rvest",
  "ggplot2", "plotly", "patchwork", "corrplot",
  "shiny", "shinydashboard", "shinyWidgets", "DT", "fresh",
  "knitr", "kableExtra",
  "tictoc", "cli", "glue"
)

install_if_needed(core_packages)

# Torch for deep learning
if(!require("torch", quietly = TRUE)) {
  cat("Installing torch (may take several minutes)...\n")
  install.packages("torch")
  torch::install_torch()
}

# Bioconductor packages
if(!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cat("\n✓ Package installation complete!\n")

# ============================================================================
# CREATE PROJECT STRUCTURE
# ============================================================================

cat("\n📁 Creating project directories...\n")

dirs <- c("data", "models", "results", "plots", "logs")
for(dir in dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("  Created:", dir, "\n")
  } else {
    cat("  Exists:", dir, "\n")
  }
}

# ============================================================================
# LOAD LIBRARIES
# ============================================================================

cat("\n📚 Loading libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rcdk)
  library(ranger)
  library(xgboost)
  library(torch)
  library(future)
  library(furrr)
  library(progressr)
  library(memuse)
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(plotly)
  library(cli)
  library(glue)
  library(tictoc)
})

cat("✓ Libraries loaded successfully!\n")

# ============================================================================
# SYSTEM DIAGNOSTICS
# ============================================================================

cat("\n🔍 System Diagnostics\n")
cat(cli::rule(), "\n")

cat("R Version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("CPU Cores:", parallel::detectCores(), "\n")

mem_info <- Sys.meminfo()
cat(sprintf("Total RAM: %.2f GB | Free RAM: %.2f GB\n",
            as.numeric(mem_info$totalram)/1e9,
            as.numeric(mem_info$freeram)/1e9))




# Check GPU
if(torch::cuda_is_available()) {
  cat(cli::col_green("✓ CUDA GPU Available!\n"))
  system("nvidia-smi", intern = TRUE)
} else {
  cat(cli::col_yellow("⚠ No GPU detected. Will use CPU.\n"))
}

cat(cli::rule(), "\n")

# ============================================================================
# CONFIGURATION OPTIONS
# ============================================================================

cat("\n⚙️  Configuration Options\n")
cat(cli::rule(), "\n")

# User configurable parameters
CONFIG <- list(
  # Data collection
  MAX_COMPOUNDS = 2000,  # Adjust based on your system
  TARGET_CHEMBL_ID = "CHEMBL3927",  # SARS-CoV-2 3CL protease
  
  # Feature calculation
  FINGERPRINT_SIZE = 1024,
  
  # Model training
  RF_NUM_TREES = 500,
  XGB_NROUNDS = 500,
  NN_EPOCHS = 100,
  NN_BATCH_SIZE = 64,
  
  # Parallel processing
  N_CORES = max(1, parallel::detectCores() - 1),
  
  # Execution control
  RUN_DATA_COLLECTION = TRUE,
  RUN_DESCRIPTOR_CALC = TRUE,
  RUN_MODEL_TRAINING = TRUE,
  LAUNCH_DASHBOARD = TRUE
)

cat("Configuration:\n")
for(name in names(CONFIG)) {
  cat(sprintf("  %s: %s\n", name, CONFIG[[name]]))
}

cat(cli::rule(), "\n")

# Confirmation
cat("\n⏱️  Estimated Total Time: 1.5 - 3 hours (depending on system)\n")
cat("💾 Required Disk Space: ~5 GB\n")
cat("🔧 System will automatically manage resources\n\n")

proceed <- readline(prompt = "Proceed with pipeline execution? (yes/no): ")

if(tolower(proceed) != "yes") {
  cat("\n❌ Pipeline execution cancelled.\n")
  quit(save = "no")
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

cat("\n🔧 Loading helper functions...\n")

source_if_exists <- function(file) {
  if(file.exists(file)) {
    source(file)
    return(TRUE)
  }
  return(FALSE)
}

# Setup parallel processing
setup_parallel <- function(cores = CONFIG$N_CORES) {
  if(Sys.getenv("RSTUDIO") == "1") {
    plan(multisession, workers = cores)
  } else {
    plan(multicore, workers = cores)
  }
  cat(glue("Parallel backend: {cores} workers\n"))
}

# Memory checker
check_memory <- function(stage = "") {
  mem <- Sys.meminfo()
  used_pct <- round((1 - as.numeric(mem$freeram) / as.numeric(mem$totalram)) * 100, 1)
  
  cat(cli::rule(glue("Memory: {stage}")), "\n")
  cat(glue("  Used: {used_pct}% | Free: {mem$freeram}\n"))
  
  if(used_pct > 85) {
    warning("Memory usage > 85%. Consider reducing dataset size.")
  }
  
  invisible(used_pct)
}

cat("✓ Helper functions loaded\n")

# ============================================================================
# PIPELINE EXECUTION LOG
# ============================================================================

log_file <- file.path("logs", paste0("pipeline_", Sys.Date(), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")

cat("\n" , cli::rule("PIPELINE EXECUTION LOG"), "\n")
cat("Started at:", as.character(Sys.time()), "\n")
cat(cli::rule(), "\n\n")

# ============================================================================
# STAGE 1: DATA COLLECTION
# ============================================================================

if(CONFIG$RUN_DATA_COLLECTION) {
  cat("\n")
  cli::cli_h1("STAGE 1: Data Collection from ChEMBL")
  
  tic("Stage 1")
  
  # Source data collection script
  if(file.exists("scripts/01_data_collection.R")) {
    source("scripts/01_data_collection.R")
  } else {
    cat("⚠ Data collection script not found. Please run QMD chunks manually.\n")
  }
  
  toc()
  check_memory("After Data Collection")
  
} else {
  cat("\n⏭️  Skipping data collection (already completed)\n")
}

# ============================================================================
# STAGE 2: MOLECULAR DESCRIPTORS
# ============================================================================

if(CONFIG$RUN_DESCRIPTOR_CALC) {
  cat("\n")
  cli::cli_h1("STAGE 2: Molecular Descriptor Calculation")
  
  tic("Stage 2")
  
  setup_parallel()
  
  if(file.exists("scripts/02_descriptors.R")) {
    source("scripts/02_descriptors.R")
  } else {
    cat("⚠ Descriptor script not found. Please run QMD chunks manually.\n")
  }
  
  toc()
  check_memory("After Descriptors")
  gc()
  
} else {
  cat("\n⏭️  Skipping descriptor calculation (already completed)\n")
}

# ============================================================================
# STAGE 3: MODEL TRAINING
# ============================================================================

if(CONFIG$RUN_MODEL_TRAINING) {
  cat("\n")
  cli::cli_h1("STAGE 3: Machine Learning Model Training")
  
  tic("Stage 3")
  
  if(file.exists("scripts/03_model_training.R")) {
    source("scripts/03_model_training.R")
  } else {
    cat("⚠ Model training script not found. Please run QMD chunks manually.\n")
  }
  
  toc()
  check_memory("After Model Training")
  gc()
  
} else {
  cat("\n⏭️  Skipping model training (already completed)\n")
}

# ============================================================================
# STAGE 4: RESULTS VISUALIZATION
# ============================================================================

cat("\n")
cli::cli_h1("STAGE 4: Generating Visualizations")

# Generate summary plots
if(file.exists("data/full_dataset_with_descriptors.csv") && 
   file.exists("results/model_comparison.csv")) {
  
  cat("Creating summary visualizations...\n")
  
  # Load data
  full_data <- read_csv("data/full_dataset_with_descriptors.csv", show_col_types = FALSE)
  comparison <- read_csv("results/model_comparison.csv", show_col_types = FALSE)
  
  # Chemical space plot
  p1 <- ggplot(full_data, aes(x = MW, y = LogP, color = bioactivity_class)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 500, linetype = "dashed", color = "red") +
    labs(title = "Chemical Space Analysis",
         x = "Molecular Weight (Da)",
         y = "LogP") +
    theme_minimal()
  
  ggsave("plots/final_chemical_space.png", p1, width = 10, height = 6, dpi = 300)
  
  # Model comparison
  p2 <- ggplot(comparison, aes(x = reorder(Model, -R_squared), y = R_squared, fill = Model)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = round(R_squared, 3)), vjust = -0.5) +
    labs(title = "Model Performance Comparison",
         x = "Model", y = "R² Score") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("plots/final_model_comparison.png", p2, width = 8, height = 6, dpi = 300)
  
  cat("✓ Plots saved to plots/ directory\n")
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cli::cli_h1("PIPELINE EXECUTION COMPLETE!")
cat(cli::rule(), "\n")

cat("\n📊 Summary of Results:\n\n")

# Check what was generated
files_generated <- c(
  "Data files" = file.exists("data/full_dataset_with_descriptors.csv"),
  "RF model" = file.exists("models/random_forest_model.rds"),
  "XGBoost model" = file.exists("models/xgboost_model.json"),
  "Neural Network" = file.exists("models/torch_nn_model.pt"),
  "Predictions" = file.exists("results/model_comparison.csv")
)

for(item in names(files_generated)) {
  if(files_generated[item]) {
    cat(cli::col_green("✓"), item, "\n")
  } else {
    cat(cli::col_red("✗"), item, "\n")
  }
}

cat("\n📁 Output Locations:\n")
cat("  Data: data/\n")
cat("  Models: models/\n")
cat("  Results: results/\n")
cat("  Plots: plots/\n")
cat("  Logs:", log_file, "\n")

cat("\n🎯 Next Steps:\n")
cat("  1. Review model_comparison.csv for performance metrics\n")
cat("  2. Launch Shiny dashboard: shiny::runApp('shiny-dashboard.R')\n")
cat("  3. Explore chemical space plots in plots/ directory\n")
cat("  4. Use trained models for new predictions\n")

# ============================================================================
# LAUNCH DASHBOARD (OPTIONAL)
# ============================================================================

if(CONFIG$LAUNCH_DASHBOARD && file.exists("shiny-dashboard.R")) {
  cat("\n")
  launch_dash <- readline(prompt = "Launch Shiny Dashboard now? (yes/no): ")
  
  if(tolower(launch_dash) == "yes") {
    cat("\n🚀 Launching interactive dashboard...\n")
    cat("Dashboard will open in your browser.\n")
    cat("Press Ctrl+C to stop the dashboard.\n\n")
    
    # Close log files before launching
    sink(type = "message")
    sink(type = "output")
    close(log_conn)
    
    shiny::runApp("shiny-dashboard.R", launch.browser = TRUE)
  }
}

# Close log file
sink(type = "message")
sink(type = "output")
close(log_conn)

cat("\n")
cat(cli::rule("Thank you for using Drug Discovery ML Pipeline!"), "\n")
cat("📧 For support, see the complete-user-guide.md\n")
cat(cli::rule(), "\n")

# End
