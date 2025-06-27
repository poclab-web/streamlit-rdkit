import pandas as pd
import numpy as np
import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

import shap

# データ読み込み
df = pd.read_csv(".data/processed/mp_logs_logp_1212_with_info.csv")

def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,))
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

os.makedirs(".models", exist_ok=True)

# 1. 融点モデル
X_mp = np.array([smiles_to_ecfp(s) for s in df["SMILES"]])
y_mp = df["MP-Measured"]
X_train, X_test, y_train, y_test = train_test_split(X_mp, y_mp, test_size=0.2, random_state=42)
mp_model = RandomForestRegressor(random_state=42)
mp_model.fit(X_train, y_train)
with open(".models/mp_rf_model.pkl", "wb") as f:
    pickle.dump(mp_model, f)
print("融点モデル保存完了")

X = np.array([smiles_to_ecfp(s) for s in df["SMILES"]])
model = pickle.load(open(".models/mp_rf_model.pkl", "rb"))

# SHAP値の計算
explainer = shap.Explainer(model, X)
shap_values = explainer(X, check_additivity=False)

# 保存
with open(".models/mp_shap_values.pkl", "wb") as f:
    pickle.dump(shap_values, f)

# 2. 溶解度モデル
y_logs = df["LogS-Measured"]
X_train, X_test, y_train, y_test = train_test_split(X_mp, y_logs, test_size=0.2, random_state=42)
logs_model = RandomForestRegressor(random_state=42)
logs_model.fit(X_train, y_train)
with open(".models/logs_rf_model.pkl", "wb") as f:
    pickle.dump(logs_model, f)
print("溶解度モデル保存完了")

# 溶解度モデルのSHAP値計算・保存
logs_explainer = shap.Explainer(logs_model, X_mp)
logs_shap_values = logs_explainer(X_mp, check_additivity=False)
with open(".models/logs_shap_values.pkl", "wb") as f:
    pickle.dump(logs_shap_values, f)
print("溶解度モデルSHAP値保存完了")

# 3. logPモデル
y_logp = df["LogP-Measured"]
X_train, X_test, y_train, y_test = train_test_split(X_mp, y_logp, test_size=0.2, random_state=42)
logp_model = RandomForestRegressor(random_state=42)
logp_model.fit(X_train, y_train)
with open(".models/logp_rf_model.pkl", "wb") as f:
    pickle.dump(logp_model, f)
print("logPモデル保存完了")

# logPモデルのSHAP値計算・保存
logp_explainer = shap.Explainer(logp_model, X_mp)
logp_shap_values = logp_explainer(X_mp, check_additivity=False)
with open(".models/logp_shap_values.pkl", "wb") as f:
    pickle.dump(logp_shap_values, f)
print("logPモデルSHAP値保存完了")