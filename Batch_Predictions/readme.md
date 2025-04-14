# HuSSPred Virtual Screening Repository

## Overview
This repository contains scripts for performing virtual screening using QSAR-ready SMILES, molecular fingerprints, and Mordred descriptors. It applies machine learning models to predict compound properties based on their chemical structures. The pipeline includes applicability domain checks, feature selection, and predictions using pre-trained models.

## Installation
### Prerequisites
Ensure you have Python (>=3.8) installed. It is recommended to use a virtual environment (e.g., `venv` or `conda`).

### Install Dependencies
#### Using pip
```bash
pip install -r requirements.txt
```

## Usage
### 1. Prepare Your Data
- Ensure your test compounds are in an Excel file with a column named **QSAR_READY_SMILES**. and in the DATA folder.

### 2. Run Virtual Screening
Execute the main script: - batch_predictions_code.ipynb

The script will:
- Read and standardize input data.
- Generate Morgan fingerprints and Mordred descriptors.
- Apply feature selection and min-max scaling.
- Load and apply pre-trained machine learning models.
- Compute applicability domain scores.
- Save results to an Excel file in the **results/** folder.

### 3. Output Files
- **2_BatchSearch_QSAR_SMILES_predictions.xlsx**: Contains binary classification results.
- **2_BatchSearch_QSAR_SMILES_predictions2.xlsx**: Contains multiclass classification results.

## Dependencies
- `numpy`
- `pandas`
- `rdkit`
- `scikit-learn`
- `joblib`
- `openpyxl`
- `mordred`

## License
This repository is open-source and freely available for research and development purposes.

## Contact
For questions, contributions, or bug reports, please open an issue on GitHub or contact the repository owner.

