# Breast Cancer Prediction with XGBoost

## Overview
This script uses XGBoost and Recursive Feature Elimination (RFE) to predict breast cancer using a genetic dataset. The model is fine-tuned with hyperparameter optimization through Grid Search and validated using cross-validation. The final predictions are saved in a submission file.

## Features
- **Feature Selection**: Uses Recursive Feature Elimination (RFE) to select the top 500 features.
- **Hyperparameter Optimization**: Optimizes XGBoost parameters via Grid Search with cross-validation.
- **Model Evaluation**: Evaluates the model using AUC-ROC scoring.
- **Prediction Output**: Generates predictions for test data and saves them to a CSV file.

## Input Data
- `train.csv`: Training dataset with the following columns:
  - `id`: Unique identifier.
  - `predicted_ancestry`: Predicted ancestry (removed in feature selection).
  - `breast_cancer`: Target variable (binary: 0 or 1).
- `test.csv`: Test dataset with the same structure as `train.csv` but without the `breast_cancer` column.

## Output
- `submission.csv`: Contains predictions for the test data:
  - `id`: Unique identifier from `test.csv`.
  - `breast_cancer`: Predicted probabilities for breast cancer.

## Steps in the Script
1. **Data Loading**:
   - Reads `train.csv` and `test.csv`.
   - Separates features (`X_train`, `X_test`) and target variable (`y_train`).

2. **Feature Selection**:
   - Applies RFE with XGBoost to select the top 500 features.

3. **Hyperparameter Optimization**:
   - Defines a hyperparameter grid.
   - Performs Grid Search with Stratified K-Fold Cross-Validation to optimize the XGBoost model.

4. **Model Training and Prediction**:
   - Fits the best model from Grid Search to the training data.
   - Generates predictions on the test data.

5. **Output**:
   - Saves the predictions to `submission.csv`.

## Hyperparameter Grid
- `learning_rate`: [0.03]
- `max_depth`: [7]
- `n_estimators`: [300]
- `subsample`: [0.8]
- `colsample_bytree`: [0.8]
- `scale_pos_weight`: Computed as the ratio of negative to positive samples.
- `reg_lambda`: [1, 10] (L2 regularization)
- `reg_alpha`: [0, 1] (L1 regularization)

## Dependencies
- `pandas`
- `numpy`
- `xgboost`
- `scikit-learn`

## Usage
Run the script with the required input files:
```bash
python XGBoost_breast_cancer_prediction.py
