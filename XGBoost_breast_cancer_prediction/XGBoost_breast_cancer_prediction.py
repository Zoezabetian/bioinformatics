import pandas as pd
from sklearn.feature_selection import RFE
from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, GridSearchCV
import numpy as np

train_data = pd.read_csv("train.csv")
test_data = pd.read_csv("test.csv")

X_train = train_data.drop(columns=["id", "predicted_ancestry", "breast_cancer"])  
y_train = train_data["breast_cancer"]

X_test = test_data.drop(columns=["id", "predicted_ancestry"])

X_test = X_test[X_train.columns]

# Recursive Feature Elimination
print("Performing Recursive Feature Elimination...")
xgb_for_rfe = XGBClassifier(random_state=42, eval_metric="auc", use_label_encoder=False)  # XGBoost as the estimator
rfe = RFE(estimator=xgb_for_rfe, n_features_to_select=500, step=50)  # select 500 features, remove 50 per iteration
X_train = rfe.fit_transform(X_train, y_train)
X_test = rfe.transform(X_test)

print(f"Selected features: {X_train.shape[1]}")

scale_pos_weight = len(y_train[y_train == 0]) / len(y_train[y_train == 1])

# hyperparameter grid
param_grid = {
    "learning_rate": [0.03],  
    "max_depth": [7],
    "n_estimators": [300],
    "subsample": [0.8],
    "colsample_bytree": [0.8],
    "scale_pos_weight": [scale_pos_weight],
    "reg_lambda": [1, 10],  # L2 regularization
    "reg_alpha": [0, 1]     # L1 regularization
}

# initialize classifier
xgb_model = XGBClassifier(random_state=42, eval_metric="auc", use_label_encoder=False)

# cross-validate
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# grid search
grid_search = GridSearchCV(
    estimator=xgb_model,
    param_grid=param_grid,
    scoring="roc_auc",
    cv=cv,
    verbose=1,
    n_jobs=-1
)
print("Starting Grid Search...")
grid_search.fit(X_train, y_train)

# best parameters and score
print("Best Parameters:", grid_search.best_params_)
print("Best Cross-Validation AUC-ROC:", grid_search.best_score_)

# use best model to predict on test data and save
best_model = grid_search.best_estimator_
y_test_pred = best_model.predict_proba(X_test)[:, 1]
submission = pd.DataFrame({
    "id": test_data["id"],  
    "breast_cancer": y_test_pred
})
submission.to_csv("submission.csv", index=False)
print("Submission file created: submission.csv")
