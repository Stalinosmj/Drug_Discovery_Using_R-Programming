# Drug Discovery Using R Programming

A comprehensive machine learning pipeline for predicting drug potency using molecular descriptors and advanced ensemble modeling techniques.

## 🎯 Project Overview

This project develops and benchmarks machine learning models to predict drug bioactivity (IC50 values) using molecular descriptors. The pipeline implements individual models (Random Forest, XGBoost, Deep Neural Networks) and advanced ensemble methods to achieve high predictive accuracy for pharmaceutical research applications.

## 📊 Key Results

- **Best Model**: Multi-Model Linear Stacking
- **Performance**: R² = 0.8974, RMSE = 0.3654, MAE = 0.2671
- **Dataset**: 2,672 compounds with 1,199 molecular features
- **Improvement**: 331% R² improvement over baseline methods

## 🔬 Dataset

- **Source**: ChEMBL Database
- **Target**: IC50 bioactivity values (converted to pChEMBL scale)
- **Total Compounds**: 2,672 unique drug-like molecules
- **Train/Test Split**: 2,140/532 compounds
- **Features**: 1,199 molecular descriptors

### Molecular Descriptors
- **Lipinski Descriptors** (9 features): MW, LogP, HBD, HBA, TPSA, nRotB, nAtoms, Aromatic Bonds, Ring Count
- **MACCS Keys** (166 features): Binary molecular fingerprints
- **ECFP4 Fingerprints** (1024 features): Circular fingerprints for structural similarity

## 🏗️ Architecture

### Individual Models
1. **Random Forest**: 1000 decision trees with bootstrap sampling
2. **XGBoost**: Gradient boosting with CPU optimization
3. **Deep Neural Network**: Multi-layer perceptron with residual connections

### Ensemble Models
1. **XGBoost + Random Forest Ensembles**:
   - Equal Weight (50-50)
   - Performance-based Weighting
   - Meta-Model Stacking
   - Bayesian Model Averaging (BMA)

2. **Multi-Model Ensembles**:
   - Weighted Average
   - Median Ensemble
   - Trimmed Mean
   - Linear Stacking (Best Performance)

## 📈 Model Performance

| Model | R² | RMSE | MAE |
|-------|-----|------|-----|
| Multi: Linear Stacking | 0.8974 | 0.3654 | 0.2671 |
| XGB+RF: Stacking | 0.8960 | 0.3678 | 0.2673 |
| XGBoost (CPU) | 0.8959 | 0.3687 | 0.2671 |
| Random Forest | 0.8743 | 0.4200 | 0.3012 |
| DNN Advanced | 0.7890 | 0.5570 | 0.4102 |

## 📊 Visualizations

The project includes 5 interactive visualizations:
1. R² Performance Ranking Chart
2. Multi-Metric Comparison
3. Speed vs Accuracy Trade-off Analysis
4. Use Case Recommendation Heatmap
5. Top 5 Models Detailed Comparison

## 🛠️ Technologies Used

- **R Programming**: Core development environment
- **Machine Learning**: `randomForest`, `xgboost`, `torch`
- **Molecular Descriptors**: `rcdk`, custom fingerprint calculation
- **Visualization**: `plotly`, `ggplot2`
- **Data Processing**: `tidyverse`, `dplyr`
- **Performance Metrics**: Custom evaluation functions

## 📋 Requirements

```r
# Required R packages
install.packages(c(
  "randomForest", "xgboost", "torch", "rcdk",
  "tidyverse", "plotly", "cli", "tictoc",
  "knitr", "readr"
))
```

## 🚀 Usage

### Basic Model Training
```r
# Load data and prepare features
source("data_preprocessing.R")
source("descriptor_calculation.R")

# Train individual models
source("train_random_forest.R")
source("train_xgboost.R")
source("train_deep_neural_network.R")

# Create ensembles
source("ensemble_models.R")

# Generate results and visualizations
source("production_pipeline_summary.R")
```

### Quick Start
```r
# Run complete pipeline
source("main_pipeline.R")
```

## 📁 Project Structure

```
Drug-Discovery-R/
├── data/
│   ├── bioactivity_clean_massive.csv
│   ├── lipinski_descriptors_massive.csv
│   └── full_dataset_massive.csv
├── models/
│   ├── rf_model_massive.rds
│   ├── xgb_model_massive.json
│   └── dnn_model_best.pt
├── results/
│   ├── final_complete_model_comparison.csv
│   ├── all_model_predictions.csv
│   └── visualizations/
├── src/
│   ├── data_preprocessing.R
│   ├── descriptor_calculation.R
│   ├── model_training.R
│   └── ensemble_methods.R
└── README.md
```

## 💡 Model Recommendations

### Use Case Scenarios

- **🚀 Fast Predictions**: XGBoost (R² = 0.8959, ~1ms inference)
- **💎 Maximum Accuracy**: Multi-Model Linear Stacking (R² = 0.8974)
- **⚖️ Balanced Performance**: XGB+RF Weighted Ensemble
- **🔬 Regulatory Applications**: Multi-Model Ensembles (most robust)

## 🔬 Research Applications

- **Drug Potency Screening**: High-throughput virtual screening
- **Lead Optimization**: Structure-activity relationship analysis  
- **ADMET Prediction**: Absorption, Distribution, Metabolism, Excretion, Toxicity
- **Chemical Space Exploration**: Novel compound discovery

## 📚 References

- Chen, H., et al. (2018). Deep learning for drug discovery. *Nature Reviews Drug Discovery*.
- Wu, Z., et al. (2018). MoleculeNet: A benchmark for molecular machine learning. *Chemical Science*.
- Mayr, A., et al. (2016). DeepTox: Toxicity prediction using deep learning. *Frontiers in Environmental Science*.

## 🤝 Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit changes (`git commit -m 'Add AmazingFeature'`)
4. Push to branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Your Name**
- GitHub: [@yourusername](https://github.com/yourusername)
- Email: your.email@example.com

## 🙏 Acknowledgments

- ChEMBL Database for providing high-quality bioactivity data
- R Community for excellent machine learning packages
- Academic institutions supporting open-source drug discovery research

---

**⭐ Star this repository if you found it helpful!**