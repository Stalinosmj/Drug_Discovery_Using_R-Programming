# Bioinformatics Drug Discovery Using Machine Learning in R

## Complete End-to-End Implementation Guide

### 📋 Table of Contents

1. [Project Overview](#project-overview)
2. [System Requirements](#system-requirements)
3. [Installation Instructions](#installation)
4. [Project Structure](#project-structure)
5. [Step-by-Step Usage Guide](#usage-guide)
6. [Features](#features)
7. [Model Details](#model-details)
8. [Troubleshooting](#troubleshooting)
9. [Performance Optimization](#optimization)
10. [References](#references)

---

## 🎯 Project Overview

This project implements a **complete drug discovery pipeline** using machine learning to predict bioactivity of chemical compounds against SARS-CoV-2 3C-like proteinase. The implementation includes:

- **Data Collection** from ChEMBL database via REST API
- **Molecular Descriptor Calculation** (Lipinski descriptors + fingerprints)
- **Multiple ML Models** (Random Forest, XGBoost, Deep Neural Networks)
- **GPU Acceleration** support via PyTorch/Torch for R
- **Parallel Processing** for faster computation
- **Progress Tracking** with time estimation
- **Resource Monitoring** to prevent system overload
- **Interactive Shiny Dashboard** for visualization and exploration

### Key Technologies

- **R 4.3+** with modern packages
- **rcdk** for cheminformatics
- **ranger** for fast Random Forest
- **xgboost** for gradient boosting
- **torch** for deep learning with GPU support
- **future/furrr** for parallel processing
- **progressr** for progress bars
- **memuse** for memory monitoring
- **Shiny** for interactive dashboards

---

## 💻 System Requirements

### Minimum Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **CPU** | 4 cores (Intel i5/Ryzen 5) | 8-12 cores (Intel i7/i9, Ryzen 7/9) |
| **RAM** | 8 GB | 16-32 GB |
| **Storage** | 10 GB free | SSD with 30+ GB |
| **GPU** | Not required | NVIDIA GPU with CUDA 11+ |
| **OS** | Windows 10/11, macOS 10.15+, Linux | Ubuntu 20.04+, macOS 12+ |

### Software Requirements

- **R** version 4.3.0 or higher
- **RStudio** 2023.06+ (recommended)
- **Java JRE** 8 or higher (for PaDEL-Descriptor/rcdk)
- **CUDA Toolkit** 11.7+ (optional, for GPU acceleration)
- **Git** (for version control)

---

## 📦 Installation Instructions

### Step 1: Install R and RStudio

1. Download and install **R** from [CRAN](https://cran.r-project.org/)
2. Download and install **RStudio** from [Posit](https://posit.co/download/rstudio-desktop/)

### Step 2: Install System Dependencies

#### On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y \
    default-jre \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev
```

#### On macOS:
```bash
brew install openjdk
brew install cairo
```

#### On Windows:
- Install Java from [Oracle](https://www.java.com/download/)
- Install Rtools from [CRAN](https://cran.r-project.org/bin/windows/Rtools/)

### Step 3: Install R Packages

Open RStudio and run:

```r
# Install all required packages
source("install_packages.R")
```

Or manually:

```r
# CRAN packages
install.packages(c(
  "tidyverse", "data.table", "rcdk", "ranger", "xgboost",
  "torch", "caret", "future", "furrr", "progressr",
  "memuse", "httr", "jsonlite", "ggplot2", "plotly",
  "shiny", "shinydashboard", "DT", "tictoc", "cli"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ChemmineR", "ChemmineOB"))

# Torch installation (with GPU support if available)
torch::install_torch()
```

### Step 4: GPU Setup (Optional but Recommended)

If you have an NVIDIA GPU:

1. Install **CUDA Toolkit** from [NVIDIA](https://developer.nvidia.com/cuda-downloads)
2. Install **cuDNN** library
3. Verify GPU availability in R:

```r
library(torch)
torch::cuda_is_available()  # Should return TRUE
torch::cuda_get_device_name(0)  # Shows your GPU name
```

### Step 5: Create Project Structure

```r
# Create necessary directories
dir.create("data", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
```

---

## 📁 Project Structure

```
drug-discovery-ml/
│
├── drug-discovery-main.qmd        # Main analysis (Part 1)
├── drug-discovery-part2.qmd       # Model training (Part 2)
├── shiny-dashboard.R              # Interactive dashboard
├── README.md                      # This file
│
├── data/
│   ├── bioactivity_raw.csv
│   ├── bioactivity_clean.csv
│   ├── full_dataset_with_descriptors.csv
│   ├── train_data.csv
│   └── test_data.csv
│
├── models/
│   ├── random_forest_model.rds
│   ├── final_rf_model_tuned.rds
│   ├── xgboost_model.json
│   └── torch_nn_model.pt
│
├── results/
│   ├── rf_test_predictions.csv
│   ├── xgb_test_predictions.csv
│   ├── nn_test_predictions.csv
│   ├── model_comparison.csv
│   ├── hyperparameter_tuning_results.csv
│   └── nn_training_history.csv
│
└── plots/
    ├── chemical_space_analysis.png
    └── bioactivity_distribution.png
```

---

## 🚀 Step-by-Step Usage Guide

### Phase 1: Data Collection and Preprocessing

#### 1.1 Open the Main QMD File

```r
# In RStudio
file.edit("drug-discovery-main.qmd")
```

#### 1.2 Run Setup and Data Collection

Execute the chunks in order:

1. **Setup and Dependencies** - Loads all packages
2. **ChEMBL API Functions** - Defines data fetching functions
3. **Fetch SARS-CoV-2 Data** - Downloads bioactivity data (takes 20-40 min)
4. **Data Preprocessing** - Cleans and prepares data

**Expected Time:** 30-60 minutes

**Output Files:**
- `data/bioactivity_raw.csv`
- `data/bioactivity_clean.csv`

#### 1.3 Calculate Molecular Descriptors

Execute descriptor calculation chunks:

1. **Lipinski Descriptors** - Calculates MW, LogP, HBD, HBA, TPSA
2. **Molecular Fingerprints** - Generates 1024-bit extended fingerprints
3. **Combine Features** - Merges all features into one dataset

**Expected Time:** 15-30 minutes (with 8 cores)

**Output Files:**
- `data/full_dataset_with_descriptors.csv`

**Progress Monitoring:**
```r
# You'll see real-time progress like:
Processing molecule 1523/3000 [=====>    ] 51% ETA: 8m
```

### Phase 2: Model Training

#### 2.1 Open Part 2 QMD File

```r
file.edit("drug-discovery-part2.qmd")
```

#### 2.2 Train Models

Execute model training chunks in sequence:

**Random Forest (Ranger)**
```r
# Expected time: 10-15 minutes
# Progress shown for each tree
# Automatic memory monitoring
```

**Hyperparameter Tuning**
```r
# Expected time: 30-45 minutes
# Tests multiple parameter combinations
# Uses 5-fold cross-validation
```

**XGBoost**
```r
# Expected time: 15-20 minutes
# Shows round-by-round RMSE
# Early stopping if no improvement
```

**Deep Neural Network (GPU)**
```r
# Expected time: 5-10 minutes (GPU) / 30-45 minutes (CPU)
# Real-time loss tracking
# Automatic GPU detection and usage
```

**Output Files:**
- `models/random_forest_model.rds`
- `models/final_rf_model_tuned.rds`
- `models/xgboost_model.json`
- `models/torch_nn_model.pt`
- All prediction CSVs in `results/`

### Phase 3: Interactive Dashboard

#### 3.1 Launch Shiny Dashboard

```r
# Run the dashboard
shiny::runApp("shiny-dashboard.R")
```

#### 3.2 Explore Dashboard Features

**Dashboard Tab:**
- Overview statistics
- Bioactivity distribution
- Model comparison

**Data Explorer:**
- Filter by bioactivity class
- Filter by pIC50 range
- View drug-like compounds only
- Interactive data table

**Chemical Space:**
- 2D scatter plots
- Customizable axes
- Lipinski rule violations

**Model Performance:**
- Predicted vs Actual plots
- Residual analysis
- Performance metrics table
- Training history (NN)

**Predictions:**
- Upload new SMILES
- Generate predictions
- Download results

**Feature Importance:**
- Top contributing features
- Compare RF and XGBoost

**System Monitor:**
- CPU usage
- Memory usage
- GPU status
- Disk usage

---

## ✨ Features

### 1. Progress Bars with Time Estimation

Every major operation shows:
- Current progress percentage
- Estimated time remaining
- Steps completed/total
- Human-readable messages

Example output:
```
[=========>              ] 45% | ETA: 12m 34s | Processing compound 1350/3000
```

### 2. Memory Management

Automatic monitoring prevents system crashes:
```
Memory Status: After descriptor calculation
  Available RAM: 8.5 GB
  Process Usage: 4.2 GB
  % Used: 47%
```

### 3. Parallel Processing

Automatically uses available cores:
```r
# Auto-detects optimal cores
setup_parallel()  
# Uses: 11 cores (out of 12, saves 1 for system)
```

### 4. GPU Acceleration

Torch neural network automatically uses GPU if available:
```
✓ CUDA GPU available!
  Device: NVIDIA GeForce RTX 3080
  Training on cuda
```

### 5. Resource Optimization

- **Saves memory**: Cleans intermediate variables
- **Efficient data structures**: Uses data.table for large datasets
- **Sparse matrices**: For fingerprint data
- **Batch processing**: Prevents memory overflow

---

## 🤖 Model Details

### Random Forest (Ranger)

**Architecture:**
- 500 trees (tunable)
- mtry = sqrt(features)
- Permutation importance

**Advantages:**
- Fast training
- Handles high dimensions well
- No hyperparameter tuning needed
- Built-in variable importance

**Performance:**
- Training time: ~10 minutes (500 trees, 8 cores)
- Typical R²: 0.65-0.75

### XGBoost

**Architecture:**
- Gradient boosting
- Learning rate: 0.1
- Max depth: 6
- Early stopping: 50 rounds

**Advantages:**
- Often best performance
- Handles missing values
- Built-in regularization
- Feature importance

**Performance:**
- Training time: ~15 minutes
- Typical R²: 0.70-0.80

### Deep Neural Network (Torch)

**Architecture:**
- Input layer: n_features
- Hidden layers: 512 → 256 → 128 → 64
- Output: 1 (pIC50 prediction)
- Dropout: 0.3, 0.3, 0.2
- Batch normalization
- ReLU activation

**Advantages:**
- GPU acceleration
- Learns complex patterns
- Flexible architecture

**Performance:**
- Training time: ~5 minutes (GPU) / ~30 minutes (CPU)
- Typical R²: 0.65-0.78

---

## 🔧 Troubleshooting

### Common Issues and Solutions

#### 1. Java Not Found (rcdk error)

**Error:**
```
Error in .jnew(...) : Java Virtual Machine not found
```

**Solution:**
```r
# Set JAVA_HOME
Sys.setenv(JAVA_HOME = "/path/to/java")
library(rJava)
.jinit()
```

#### 2. Out of Memory

**Error:**
```
Error: cannot allocate vector of size X GB
```

**Solution:**
```r
# Reduce dataset size
bioactivity_clean <- bioactivity_clean %>% slice_sample(n = 2000)

# Or increase memory limit (Windows)
memory.limit(size = 16000)  # 16 GB
```

#### 3. CUDA Not Available

**Error:**
```
⚠ No CUDA GPU detected. Using CPU.
```

**Solution:**
```r
# Reinstall torch with CUDA support
torch::install_torch(reinstall = TRUE)

# Check CUDA version compatibility
system("nvidia-smi")
```

#### 4. Slow Descriptor Calculation

**Solution:**
```r
# Reduce fingerprint size
molecular_fingerprints <- calculate_fingerprints(
  smiles_vector,
  fp_type = "extended",
  size = 512  # Instead of 1024
)

# Use fewer cores if system is overloaded
plan(multisession, workers = 4)
```

#### 5. ChEMBL API Timeout

**Error:**
```
API error: Timeout was reached
```

**Solution:**
```r
# Increase timeout
response <- GET(
  base_url,
  query = list(...),
  timeout(60)  # Increase from 30
)

# Reduce batch size
limit <- 500  # Instead of 1000
```

---

## ⚡ Performance Optimization

### Tips for Faster Execution

#### 1. Use SSD Storage
```r
# Store data on SSD, not HDD
setwd("/path/to/ssd/drug-discovery-ml")
```

#### 2. Optimize Parallel Workers
```r
# Leave 1-2 cores for system
cores <- parallel::detectCores() - 2
plan(multisession, workers = cores)
```

#### 3. Use ranger Instead of randomForest
```r
# ranger is 10x faster
library(ranger)  # ✓ Fast
# library(randomForest)  # ✗ Slow
```

#### 4. Enable GPU for Neural Network
```r
# 10-20x speedup
device <- torch_device("cuda")
```

#### 5. Reduce Dataset for Testing
```r
# Test with smaller dataset first
sample_data <- full_dataset %>% slice_sample(n = 500)
```

#### 6. Save Intermediate Results
```r
# Avoid re-computation
write_csv(descriptors, "data/descriptors_cache.csv")
```

---

## 📊 Expected Runtime Summary

| Task | Laptop (8GB, 4 cores) | Workstation (32GB, 12 cores, GPU) |
|------|----------------------|-----------------------------------|
| **Data Collection** | 40 min | 30 min |
| **Descriptor Calculation** | 60 min | 15 min |
| **Random Forest Training** | 25 min | 8 min |
| **Hyperparameter Tuning** | 90 min | 35 min |
| **XGBoost Training** | 30 min | 12 min |
| **Neural Network (CPU)** | 45 min | - |
| **Neural Network (GPU)** | - | 5 min |
| **Total (without tuning)** | 3.5 hours | 1.2 hours |
| **Total (with tuning)** | 5 hours | 2 hours |

---

## 📚 References

1. **ChEMBL Database**  
   Gaulton et al. (2024). "ChEMBL Database in 2023." Nucleic Acids Research.  
   https://www.ebi.ac.uk/chembl/

2. **Lipinski's Rule of Five**  
   Lipinski et al. (2001). "Experimental and computational approaches to estimate solubility and permeability."  
   Advanced Drug Delivery Reviews.

3. **PaDEL-Descriptor**  
   Yap (2011). "PaDEL-descriptor: An open source software to calculate molecular descriptors and fingerprints."  
   Journal of Computational Chemistry.

4. **Random Forest (Ranger)**  
   Wright & Ziegler (2017). "ranger: A Fast Implementation of Random Forests."  
   https://github.com/imbs-hl/ranger

5. **XGBoost**  
   Chen & Guestrin (2016). "XGBoost: A Scalable Tree Boosting System."  
   KDD '16 Proceedings.

6. **Torch for R**  
   Falbel & Luraschi (2021). "torch: Tensors and Neural Networks with GPU Acceleration."  
   https://torch.mlverse.org/

---

## 👨‍💻 Author and License

**Author:** Drug Discovery ML Team  
**License:** MIT  
**Version:** 1.0.0  
**Last Updated:** October 2025

---

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## ❓ Support

For questions or issues:
- Open a GitHub issue
- Check troubleshooting section
- Review R package documentation

---

**Happy Drug Discovery! 🧬💊**
