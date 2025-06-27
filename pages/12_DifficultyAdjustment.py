import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt  # 追加
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import Draw 
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

import shap
import pickle

from logic.FingerPrint import morgan_fingerprint


# アプリの定義

df = pd.read_csv(".models/mp_logs_logp_1212_with_info.csv")

def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((1,))
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

@st.cache_resource
def load_rf_model(model_path):
    with open(model_path, "rb") as f:
        return pickle.load(f)

@st.cache_resource
def load_shap_values(path):
    with open(path, "rb") as f:
        return pickle.load(f)

@st.cache_resource
def load_gse_logs():
    return pd.read_csv(".models/logs_gse.csv")["LogS-GSE"]

def difficulty_adjustment_display():
    """
    融点、溶解度、logPデータの確認
    """
    st.title("融点、溶解度、logPデータの確認")
    st.write("融点、溶解度、logPのデータを確認します。")

    # データフレームの表示
    st.subheader("データフレーム")
    st.dataframe(df)

    # 散布図行列の表示
    st.subheader("散布図行列")
    pairplot = sns.pairplot(df)
    st.pyplot(pairplot.figure)

    # インデックスで分子を選択
    st.subheader("分子の選択")
    selected_idx = st.number_input("分子のindexを選択してください", min_value=0, max_value=len(df)-1, value=0, step=1)
    selected_row = df.iloc[int(selected_idx)]
    selected_smiles = selected_row["SMILES"]


    # 分子の情報を表示
    st.subheader("分子の情報")
    mol = Chem.MolFromSmiles(selected_smiles)
    if mol:
        st.write(f"SMILES: {selected_smiles}")
        st.write(f"分子量: {Chem.Descriptors.MolWt(mol):.2f}")
        st.write(f"融点（実測）: {selected_row['MP-Measured']}")
        st.write(f"溶解度（実測）: {selected_row['LogS-Measured']}")
        st.write(f"logP（実測）: {selected_row['LogP-Measured']}")
    else:
        st.error("無効なSMILESです。")
        return

    # fingerprintのbit画像表示
    st.subheader("分子のFingerprint表示")

    img, bitI_morgan, mol = morgan_fingerprint(selected_smiles, radius=2, nBits=2048)

    st.text("1となっているbitの合計数(radius=2, nBits=2048)")
    st.code(len(bitI_morgan.keys()))

    st.text("1となっている部分構造を表示")
    st.image(img)

    st.text("分子の構造を表示")

    for atom in mol.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

    img2 = Draw.MolToImage(mol)
    st.image(img2)

    st.text("それぞれのbitに入っている情報を表示 IDを表示させた後 （中心原子の番号, radius) を表示させる")
    for bit, value in bitI_morgan.items():
        st.text(bit)
        st.text(value)



def melting_point_prediction():
    """
    融点予測のページ
    """
    st.subheader("ランダムフォレスト回帰による融点予測（ECFP特徴量）")

    # ECFP特徴量の作成
    ecfp_features = np.array([smiles_to_ecfp(s) for s in df["SMILES"]])
    X = ecfp_features
    y = df["MP-Measured"]

    # データ分割
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # モデル読み込み
    model = load_rf_model(".models/mp_rf_model.pkl")

    # 予測
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # 評価
    mse_train = mean_squared_error(y_train, y_train_pred)
    rmse_train = np.sqrt(mse_train)
    mae_train = np.mean(np.abs(y_train - y_train_pred))
    r2_train = r2_score(y_train, y_train_pred)

    mse_test = mean_squared_error(y_test, y_test_pred)
    rmse_test = np.sqrt(mse_test)
    mae_test = np.mean(np.abs(y_test - y_test_pred))
    r2_test = r2_score(y_test, y_test_pred)

    # 評価指標の表示
    st.subheader("📈 評価指標")
    st.markdown(f"**訓練データ** - MSE: {mse_train:.3f}, RMSE: {rmse_train:.3f}, MAE: {mae_train:.3f}, R²: {r2_train:.3f}")
    st.markdown(f"**テストデータ** - MSE: {mse_test:.3f}, RMSE: {rmse_test:.3f}, MAE: {mae_test:.3f}, R²: {r2_test:.3f}")

    # yyプロット（訓練＋テスト）
    st.subheader("🔍 yyプロット（実測値 vs 予測値）")
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_train_pred, label="train", color="blue", alpha=0.6)
    ax.scatter(y_test, y_test_pred, label="test", color="red", alpha=0.6)
    ax.plot([min(y), max(y)], [min(y), max(y)], 'k--', lw=2)
    ax.set_xlabel("actual")
    ax.set_ylabel("predict")
    ax.set_title("yyplot")
    ax.legend()
    st.pyplot(fig)

    # SHAP値の読み込み
    shap_values = load_shap_values(".models/mp_shap_values.pkl")

    # データ選択
    st.subheader("解析したいデータポイントを選んでください")
    selected_idx = st.number_input(
        "分子のindexを選択してください", 
        min_value=0, 
        max_value=len(df)-1, 
        value=0, 
        step=1,
        key="mp_pred_index"  # ← ユニークなkeyを追加
    )
    
    with st.expander("選択した分子のbitの詳細"):
        selected_row = df.iloc[int(selected_idx)]
        selected_smiles = selected_row["SMILES"]

        img, bitI_morgan, mol = morgan_fingerprint(selected_smiles, radius=2, nBits=2048)

        st.text("1となっているbitの合計数(radius=2, nBits=2048)")
        st.code(len(bitI_morgan.keys()))

        st.text("1となっている部分構造を表示")
        st.image(img)

    # 特徴量表示
    st.write("### 入力特徴量と値")

    # 選択したサンプルの特徴量を取得
    selected_sample = smiles_to_ecfp(selected_smiles)
    # 予測と実測
    pred = model.predict([selected_sample])[0]
    true = selected_row["MP-Measured"]
    st.markdown(f"**実測値（融点）**: {true:.3f}")
    st.markdown(f"**予測値**: {pred:.3f}")

    # force風 summary plot
    st.subheader("特徴量ごとのSHAP値（予測への貢献度）")
    st.write("SHAP値は、各特徴量が予測にどのように寄与しているかを示します。")
    st.write("SHAP値が正の値であれば、その特徴量は予測値を上げる方向に寄与し、負の値であれば下げる方向に寄与します。")
    st.write("SHAP値の絶対値が大きいほど、その特徴量の予測への影響が大きいことを示します。")    
    fig, ax = plt.subplots()
    shap.plots.bar(shap_values[int(selected_idx)], show=False)
    st.pyplot(fig)

    # SHAP summary plot
    st.subheader("① Summary Plotで全体の傾向を見る")
    fig_summary, ax_summary = plt.subplots()
    shap.plots.beeswarm(shap_values, show=False)
    st.pyplot(fig_summary)



def solubility_prediction_display():
    """
    溶解度予測のページ
    """
    st.title("溶解度予測")
    ecfp_features = np.array([smiles_to_ecfp(s) for s in df["SMILES"]])
    X = ecfp_features
    y = df["LogS-Measured"]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # モデル読み込み
    model = load_rf_model(".models/logs_rf_model.pkl")

    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # 評価
    mse_train = mean_squared_error(y_train, y_train_pred)
    rmse_train = np.sqrt(mse_train)
    mae_train = np.mean(np.abs(y_train - y_train_pred))
    r2_train = r2_score(y_train, y_train_pred)

    mse_test = mean_squared_error(y_test, y_test_pred)
    rmse_test = np.sqrt(mse_test)
    mae_test = np.mean(np.abs(y_test - y_test_pred))
    r2_test = r2_score(y_test, y_test_pred)

    # 評価指標の表示
    st.subheader("📈 評価指標")
    st.markdown(f"**訓練データ** - MSE: {mse_train:.3f}, RMSE: {rmse_train:.3f}, MAE: {mae_train:.3f}, R²: {r2_train:.3f}")
    st.markdown(f"**テストデータ** - MSE: {mse_test:.3f}, RMSE: {rmse_test:.3f}, MAE: {mae_test:.3f}, R²: {r2_test:.3f}")

    # yyプロット（訓練＋テスト）
    st.subheader("🔍 yyプロット（実測値 vs 予測値）")
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_train_pred, label="train", color="blue", alpha=0.6)
    ax.scatter(y_test, y_test_pred, label="test", color="red", alpha=0.6)
    ax.plot([min(y), max(y)], [min(y), max(y)], 'k--', lw=2)
    ax.set_xlabel("actual")
    ax.set_ylabel("predict")
    ax.set_title("yyplot")
    ax.legend()
    st.pyplot(fig)

    # SHAP値の読み込みと表示
    shap_values = load_shap_values(".models/logs_shap_values.pkl")
    st.subheader("SHAP値による特徴量重要度（溶解度モデル）")
    idx = st.number_input("SHAP値を見たい分子のindexを選択してください", min_value=0, max_value=len(df)-1, value=0, step=1, key="logs_shap_index")
    fig_shap, ax_shap = plt.subplots()

    with st.expander("選択した分子のbitの詳細"):
        selected_row = df.iloc[int(idx)]
        selected_smiles = selected_row["SMILES"]

        img, bitI_morgan, mol = morgan_fingerprint(selected_smiles, radius=2, nBits=2048)

        st.text("1となっているbitの合計数(radius=2, nBits=2048)")
        st.code(len(bitI_morgan.keys()))

        st.text("1となっている部分構造を表示")
        st.image(img)


    shap.plots.bar(shap_values[int(idx)], show=False)
    st.pyplot(fig_shap)
    st.subheader("Summary Plot（全体傾向）")
    fig_summary, ax_summary = plt.subplots()
    shap.plots.beeswarm(shap_values, show=False)
    st.pyplot(fig_summary)

    # General Solubility Equation (GSE) による溶解度の計算
    st.title("GSEによる溶解度の計算")
    st.write("GSEを使用して、融点とlogPから溶解度を計算します。")
    def calculate_gse_logS(melting_point, logP):
        """
        General Solubility Equation (GSE) による溶解度(logS)の計算
        :param melting_point: 融点（℃）
        :param logP: logP値
        :return: logS値
        """
        return 0.5 - 0.01 * (melting_point - 25) - logP 
    
        # GSEによるlogS計算
    df["LogS-GSE"] = df.apply(lambda row: calculate_gse_logS(row["MP-Measured"], row["LogP-Measured"]), axis=1)
    y_true = df["LogS-Measured"]
    y_pred = df["LogS-GSE"]

    # 評価指標
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = np.mean(np.abs(y_true - y_pred))
    r2 = r2_score(y_true, y_pred)

    st.subheader("📈 GSE計算値の評価指標")
    st.markdown(f"**MSE:** {mse:.3f}, **RMSE:** {rmse:.3f}, **MAE:** {mae:.3f}, **R²:** {r2:.3f}")

    # yyプロット
    st.subheader("🔍 GSE計算値のyyプロット（実測値 vs GSE計算値）")
    fig, ax = plt.subplots()
    ax.scatter(y_true, y_pred, color="purple", alpha=0.6)
    ax.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'k--', lw=2)
    ax.set_xlabel("actual")
    ax.set_ylabel("GSE predict")
    ax.set_title("yyplot (GSE logS)")
    st.pyplot(fig)
    

    
def logp_prediction_display():
    """
    logP予測のページ
    """
    ecfp_features = np.array([smiles_to_ecfp(s) for s in df["SMILES"]])
    X = ecfp_features
    y = df["LogP-Measured"]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # モデル読み込み
    model = load_rf_model(".models/logp_rf_model.pkl")

    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # 評価
    mse_train = mean_squared_error(y_train, y_train_pred)
    rmse_train = np.sqrt(mse_train)
    mae_train = np.mean(np.abs(y_train - y_train_pred))
    r2_train = r2_score(y_train, y_train_pred)

    mse_test = mean_squared_error(y_test, y_test_pred)
    rmse_test = np.sqrt(mse_test)
    mae_test = np.mean(np.abs(y_test - y_test_pred))
    r2_test = r2_score(y_test, y_test_pred)

    # ページタイトル
    st.title("ランダムフォレストによるlogP予測")

    # 評価指標の表示
    st.subheader("📈 評価指標")
    st.markdown(f"**訓練データ** - MSE: {mse_train:.3f}, RMSE: {rmse_train:.3f}, MAE: {mae_train:.3f}, R²: {r2_train:.3f}")
    st.markdown(f"**テストデータ** - MSE: {mse_test:.3f}, RMSE: {rmse_test:.3f}, MAE: {mae_test:.3f}, R²: {r2_test:.3f}")

    # yyプロット（訓練＋テスト）
    st.subheader("🔍 yyプロット（実測値 vs 予測値）")
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_train_pred, label="train", color="blue", alpha=0.6)
    ax.scatter(y_test, y_test_pred, label="test", color="red", alpha=0.6)
    ax.plot([min(y), max(y)], [min(y), max(y)], 'k--', lw=2)
    ax.set_xlabel("actual")
    ax.set_ylabel("predict")
    ax.set_title("yyplot")
    ax.legend()
    st.pyplot(fig)

    # SHAP値の読み込みと表示
    shap_values = load_shap_values(".models/logp_shap_values.pkl")
    st.subheader("SHAP値による特徴量重要度（logPモデル）")
    idx = st.number_input("SHAP値を見たい分子のindexを選択してください", min_value=0, max_value=len(df)-1, value=0, step=1, key="logp_shap_index")

    with st.expander("選択した分子のbitの詳細"):
        selected_row = df.iloc[int(idx)]
        selected_smiles = selected_row["SMILES"]

        img, bitI_morgan, mol = morgan_fingerprint(selected_smiles, radius=2, nBits=2048)

        st.text("1となっているbitの合計数(radius=2, nBits=2048)")
        st.code(len(bitI_morgan.keys()))

        st.text("1となっている部分構造を表示")
        st.image(img)

    fig_shap, ax_shap = plt.subplots()
    shap.plots.bar(shap_values[int(idx)], show=False)
    st.pyplot(fig_shap)
    st.subheader("Summary Plot（全体傾向）")
    fig_summary, ax_summary = plt.subplots()
    shap.plots.beeswarm(shap_values, show=False)
    st.pyplot(fig_summary)

    # RDKitでlogPを計算
    rdkit_logp = df["SMILES"].apply(lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s) if Chem.MolFromSmiles(s) else Chem.MolFromSmiles('')))
    y_true = df["LogP-Measured"]
    y_pred = rdkit_logp

    # 評価指標
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = np.mean(np.abs(y_true - y_pred))
    r2 = r2_score(y_true, y_pred)

    st.subheader("📈 RDKit logP計算値の評価指標")
    st.markdown(f"**MSE:** {mse:.3f}, **RMSE:** {rmse:.3f}, **MAE:** {mae:.3f}, **R²:** {r2:.3f}")

    # yyプロット
    st.subheader("🔍 RDKit logP計算値のyyプロット（実測値 vs RDKit計算値）")
    fig2, ax2 = plt.subplots()
    ax2.scatter(y_true, y_pred, color="green", alpha=0.6)
    ax2.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], 'k--', lw=2)
    ax2.set_xlabel("actual")
    ax2.set_ylabel("RDKit predict")
    ax2.set_title("yyplot (RDKit logP)")
    st.pyplot(fig2)

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "DifficultyAdjustment"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()