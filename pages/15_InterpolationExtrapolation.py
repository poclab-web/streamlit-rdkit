import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Draw
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import os
from io import BytesIO
from PIL import Image

# ========== 1. 記述子の計算 ==========
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [np.nan] * 9

    try:
        mol = Chem.AddHs(mol)
        ComputeGasteigerCharges(mol)
        charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms()]
        total_charge = sum(charges)
        max_charge = max(charges)
    except:
        total_charge = np.nan
        max_charge = np.nan

    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        Lipinski.NumRotatableBonds(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        rdMolDescriptors.CalcNumRings(mol),
        total_charge,
        max_charge
    ]

descriptor_names = [
    'MolWt', 'LogP', 'TPSA', 'RotBonds',
    'HDonors', 'HAcceptors', 'RingCount',
    'TotalCharge', 'MaxCharge'
]

# ========== 2. 分子構造描画 ==========
def draw_molecule(smiles, size=(300, 300)):
    """SMILESから分子構造を描画"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # 分子構造を描画
    img = Draw.MolToImage(mol, size=size)
    return img

def display_molecule_with_info(smiles, scaffold=None, actual=None, predicted=None, residual=None):
    """分子構造と情報を表示"""
    img = draw_molecule(smiles)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        if img is not None:
            st.image(img, caption="分子構造", use_container_width=True)
        else:
            st.error("分子構造を描画できませんでした")
    
    with col2:
        st.write(f"**SMILES:** {smiles}")
        if scaffold is not None:
            st.write(f"**骨格:** {scaffold}")
        if actual is not None:
            st.write(f"**実測値:** {actual:.3f}")
        if predicted is not None:
            st.write(f"**予測値:** {predicted:.3f}")
        if residual is not None:
            st.write(f"**残差:** {residual:.3f}")

# ========== 3. データ読み込みと記述子追加 ==========
@st.cache_data
def load_and_process_data():
    # Excel ファイルを読み込む
    data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'cbs.xlsx')
    df = pd.read_excel(data_path)
    
    # 記述子を計算
    desc = df['SMILES'].apply(compute_descriptors)
    desc_df = pd.DataFrame(desc.tolist(), columns=descriptor_names)
    df = pd.concat([df, desc_df], axis=1)
    
    # 欠損値を削除
    df_clean = df.dropna(subset=descriptor_names + ['ΔΔG.expt.']).copy()
    
    return df_clean

# ========== 4. モデル構築と評価 ==========
def build_model(df):
    X = df[descriptor_names]
    y = df['ΔΔG.expt.']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    rf_model = RandomForestRegressor(random_state=42)
    rf_model.fit(X_train, y_train)
    y_pred_train = rf_model.predict(X_train)
    y_pred_test = rf_model.predict(X_test)
    
    return rf_model, X_train, X_test, y_train, y_test, y_pred_train, y_pred_test

# ========== 5. 可視化 ==========
def plot_predictions(y_train, y_pred_train, y_test, y_pred_test):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 実測値 vs 予測値
    ax1.scatter(y_train, y_pred_train, label='Training', alpha=0.6)
    ax1.scatter(y_test, y_pred_test, label='Test', alpha=0.8)
    min_val = min(y_train.min(), y_test.min())
    max_val = max(y_train.max(), y_test.max())
    ax1.plot([min_val, max_val], [min_val, max_val], 'k--')
    ax1.set_xlabel('Actual ΔΔG.expt.')
    ax1.set_ylabel('Predicted ΔΔG.expt.')
    ax1.set_title('Actual vs Predicted')
    ax1.legend()
    ax1.grid(True)
    
    # 残差プロット
    residuals = y_test - y_pred_test
    ax2.scatter(y_pred_test, residuals)
    ax2.axhline(0, color='black', linestyle='--')
    ax2.set_xlabel('Predicted ΔΔG.expt.')
    ax2.set_ylabel('Residuals')
    ax2.set_title('Residual Plot')
    ax2.grid(True)
    
    plt.tight_layout()
    return fig

# ========== Streamlit アプリのメイン部分 ==========
def random_split_analysis():
    st.title("機械学習を用いた不斉触媒反応の解析")
    st.markdown("化学的記述子を用いた機械学習モデルによる補間・外挿解析")
    
    # データの読み込み
    with st.spinner("データを読み込み中..."):
        df = load_and_process_data()
    
    st.success(f"データ読み込み完了: {len(df)} サンプル")
    
    # データの概要を表示
    st.subheader("データの概要")
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("サンプル数", len(df))
        st.metric("記述子数", len(descriptor_names))
    
    with col2:
        st.metric("ΔΔG.expt. 平均", f"{df['ΔΔG.expt.'].mean():.3f}")
        st.metric("ΔΔG.expt. 標準偏差", f"{df['ΔΔG.expt.'].std():.3f}")
    
    # データの一部を表示
    st.subheader("データサンプル")
    st.dataframe(df[['SMILES', 'ΔΔG.expt.'] + descriptor_names].head(10))
    
    # モデル構築
    st.subheader("モデル構築と評価")
    
    if st.button("モデルを構築"):
        with st.spinner("モデルを構築中..."):
            rf_model, X_train, X_test, y_train, y_test, y_pred_train, y_pred_test = build_model(df)
        
        # 評価指標の表示
        st.subheader("評価指標")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("テストデータ MSE", f"{mean_squared_error(y_test, y_pred_test):.3f}")
            st.metric("テストデータ R²", f"{r2_score(y_test, y_pred_test):.3f}")
        
        with col2:
            st.metric("トレーニングデータ MSE", f"{mean_squared_error(y_train, y_pred_train):.3f}")
            st.metric("トレーニングデータ R²", f"{r2_score(y_train, y_pred_train):.3f}")
        
        # モデル情報
        st.subheader("モデル情報")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**使用したモデル**")
            st.write("Random Forest Regressor")
            st.write("**パラメータ**")
            st.write(f"- n_estimators: {rf_model.n_estimators}")
            st.write(f"- max_depth: {rf_model.max_depth}")
            st.write(f"- min_samples_split: {rf_model.min_samples_split}")
            st.write(f"- min_samples_leaf: {rf_model.min_samples_leaf}")
            st.write(f"- random_state: {rf_model.random_state}")
        
        with col2:
            st.write("**データ分割情報**")
            st.write("ランダム分割")
            st.write(f"- 学習データ: {len(X_train)} サンプル ({len(X_train)/(len(X_train)+len(X_test))*100:.1f}%)")
            st.write(f"- テストデータ: {len(X_test)} サンプル ({len(X_test)/(len(X_train)+len(X_test))*100:.1f}%)")
            st.write(f"- 使用記述子数: {len(descriptor_names)}")
            st.write(f"- random_state: 42")
        
        # 可視化
        st.subheader("可視化")
        fig = plot_predictions(y_train, y_pred_train, y_test, y_pred_test)
        st.pyplot(fig)
        
        # 特徴量重要度
        st.subheader("特徴量重要度")
        feature_importance = pd.DataFrame({
            'Feature': descriptor_names,
            'Importance': rf_model.feature_importances_
        }).sort_values('Importance', ascending=False)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(feature_importance['Feature'], feature_importance['Importance'])
        ax.set_xlabel('Importance')
        ax.set_title('Feature Importance')
        plt.tight_layout()
        st.pyplot(fig)
        
        # 詳細な結果
        st.subheader("詳細な結果")
        results_df = pd.DataFrame({
            'Index': X_test.index,
            'SMILES': df.loc[X_test.index, 'SMILES'],
            'Actual': y_test,
            'Predicted': y_pred_test,
            'Residual': y_test - y_pred_test,
            'Abs_Residual': np.abs(y_test - y_pred_test)
        })
        st.dataframe(results_df)
        
        # 最も予測困難な構造の表示
        st.subheader("最も予測困難な構造")
        st.write("**誤差の大きい順 (上位5件)**")
        error_analysis = results_df.nlargest(5, 'Abs_Residual').copy()
        
        for i, (idx, row) in enumerate(error_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                None,
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()
        
        st.write("**最も予測精度の高い構造 (誤差の小さい順)**")
        accuracy_analysis = results_df.nsmallest(5, 'Abs_Residual').copy()
        
        for i, (idx, row) in enumerate(accuracy_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                None,
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()

def high_selectivity_prediction():
    st.title("機械学習を用いた不斉触媒反応の解析(高い選択性を予測)")
    st.markdown("化学的記述子を用いた機械学習モデルによる補間・外挿解析（ΔΔG.expt.の昇順に並べて分割：下位80%学習）")
    
    # データの読み込み
    with st.spinner("データを読み込み中..."):
        df = load_and_process_data()
    
    st.success(f"データ読み込み完了: {len(df)} サンプル")
    
    # データの概要を表示
    st.subheader("データの概要")
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("サンプル数", len(df))
        st.metric("記述子数", len(descriptor_names))
    
    with col2:
        st.metric("ΔΔG.expt. 平均", f"{df['ΔΔG.expt.'].mean():.3f}")
        st.metric("ΔΔG.expt. 標準偏差", f"{df['ΔΔG.expt.'].std():.3f}")
    
    # データの一部を表示
    st.subheader("データサンプル")
    st.dataframe(df[['SMILES', 'ΔΔG.expt.'] + descriptor_names].head(10))
    
    # モデル構築
    st.subheader("モデル構築と評価")
    
    # 分割方法の説明
    st.info("⚠️ 分割方法: ΔΔG.expt.の昇順に並べて、下位80%を学習データ、上位20%を予測データとして使用")
    
    if st.button("モデルを構築", key="build_model_2"):
        with st.spinner("モデルを構築中..."):
            # ΔΔG.expt.の昇順に並べて分割（下位80%学習）
            df_sorted = df.sort_values(by='ΔΔG.expt.').reset_index(drop=True)
            cutoff = int(len(df_sorted) * 0.8)
            df_train = df_sorted.iloc[:cutoff]
            df_test = df_sorted.iloc[cutoff:]
            
            X_train = df_train[descriptor_names]
            y_train = df_train['ΔΔG.expt.']
            X_test = df_test[descriptor_names]
            y_test = df_test['ΔΔG.expt.']
            
            # モデル構築
            rf_model = RandomForestRegressor(random_state=42)
            rf_model.fit(X_train, y_train)
            y_pred_train = rf_model.predict(X_train)
            y_pred_test = rf_model.predict(X_test)
        
        # 評価指標の表示
        st.subheader("評価指標")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("テストデータ MSE", f"{mean_squared_error(y_test, y_pred_test):.3f}")
            st.metric("テストデータ R²", f"{r2_score(y_test, y_pred_test):.3f}")
        
        with col2:
            st.metric("トレーニングデータ MSE", f"{mean_squared_error(y_train, y_pred_train):.3f}")
            st.metric("トレーニングデータ R²", f"{r2_score(y_train, y_pred_train):.3f}")
        
        # 分割データの統計情報
        st.subheader("分割データの統計情報")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**テストデータ (上位20%)**")
            st.metric("サンプル数", len(df_test))
            st.metric("ΔΔG.expt. 範囲", f"{y_test.min():.3f} ~ {y_test.max():.3f}")
            st.metric("ΔΔG.expt. 平均", f"{y_test.mean():.3f}")
        
        with col2:
            st.write("**学習データ (下位80%)**")
            st.metric("サンプル数", len(df_train))
            st.metric("ΔΔG.expt. 範囲", f"{y_train.min():.3f} ~ {y_train.max():.3f}")
            st.metric("ΔΔG.expt. 平均", f"{y_train.mean():.3f}")
        
        # モデル情報
        st.subheader("モデル情報")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**使用したモデル**")
            st.write("Random Forest Regressor")
            st.write("**パラメータ**")
            st.write(f"- n_estimators: {rf_model.n_estimators}")
            st.write(f"- max_depth: {rf_model.max_depth}")
            st.write(f"- min_samples_split: {rf_model.min_samples_split}")
            st.write(f"- min_samples_leaf: {rf_model.min_samples_leaf}")
            st.write(f"- random_state: {rf_model.random_state}")
        
        with col2:
            st.write("**データ分割情報**")
            st.write("ΔΔG.expt.値による分割")
            st.write(f"- 学習データ: {len(X_train)} サンプル (下位80%)")
            st.write(f"- テストデータ: {len(X_test)} サンプル (上位20%)")
            st.write(f"- 使用記述子数: {len(descriptor_names)}")
            st.write(f"- 分割基準: ΔΔG.expt.昇順")
        
        # 可視化
        st.subheader("可視化")
        fig = plot_predictions(y_train, y_pred_train, y_test, y_pred_test)
        st.pyplot(fig)
        
        # 特徴量重要度
        st.subheader("特徴量重要度")
        feature_importance = pd.DataFrame({
            'Feature': descriptor_names,
            'Importance': rf_model.feature_importances_
        }).sort_values('Importance', ascending=False)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(feature_importance['Feature'], feature_importance['Importance'])
        ax.set_xlabel('Importance')
        ax.set_title('Feature Importance')
        plt.tight_layout()
        st.pyplot(fig)
        
        # 詳細な結果
        st.subheader("詳細な結果")
        results_df = pd.DataFrame({
            'Index': X_test.index,
            'SMILES': df.loc[X_test.index, 'SMILES'],
            'Actual': y_test,
            'Predicted': y_pred_test,
            'Residual': y_test - y_pred_test,
            'Abs_Residual': np.abs(y_test - y_pred_test)
        })
        st.dataframe(results_df)
        
        # 外挿性能の評価
        st.subheader("外挿性能の評価")
        st.write("高い選択性（高いΔΔG.expt.値）を持つサンプルに対する予測性能:")
        
        # 予測誤差の分析
        mae = np.mean(np.abs(y_test - y_pred_test))
        st.metric("平均絶対誤差 (MAE)", f"{mae:.3f}")
        
        # 予測精度の分布
        accuracy_within_1 = np.mean(np.abs(y_test - y_pred_test) < 1.0) * 100
        accuracy_within_05 = np.mean(np.abs(y_test - y_pred_test) < 0.5) * 100
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("誤差±1.0以内の予測精度", f"{accuracy_within_1:.1f}%")
        with col2:
            st.metric("誤差±0.5以内の予測精度", f"{accuracy_within_05:.1f}%")
        
        # 最も予測困難な構造の表示
        st.subheader("最も予測困難な構造")
        st.write("**誤差の大きい順 (上位5件)**")
        error_analysis = results_df.nlargest(5, 'Abs_Residual').copy()
        
        for i, (idx, row) in enumerate(error_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                None,
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()
        
        st.write("**最も予測精度の高い構造 (誤差の小さい順)**")
        accuracy_analysis = results_df.nsmallest(5, 'Abs_Residual').copy()
        
        for i, (idx, row) in enumerate(accuracy_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                None,
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()

def scaffold_based_split_analysis():
    st.title("機械学習を用いた不斉触媒反応の解析(骨格分割による外挿評価)")
    st.markdown("化学的記述子を用いた機械学習モデルによる補間・外挿解析（ユニーク骨格での評価）")
    
    # データの読み込み
    with st.spinner("データを読み込み中..."):
        df = load_and_process_data()
    
    st.success(f"データ読み込み完了: {len(df)} サンプル")
    
    # データの概要を表示
    st.subheader("データの概要")
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("サンプル数", len(df))
        st.metric("記述子数", len(descriptor_names))
    
    with col2:
        st.metric("ΔΔG.expt. 平均", f"{df['ΔΔG.expt.'].mean():.3f}")
        st.metric("ΔΔG.expt. 標準偏差", f"{df['ΔΔG.expt.'].std():.3f}")
    
    # データの一部を表示
    st.subheader("データサンプル")
    st.dataframe(df[['SMILES', 'ΔΔG.expt.'] + descriptor_names].head(10))
    
    # モデル構築
    st.subheader("モデル構築と評価")
    
    # 分割方法の説明
    st.info("⚠️ 分割方法: ユニーク骨格（1件のみ）を予測データ、その他を学習データとして使用")
    
    if st.button("モデルを構築", key="build_model_3"):
        with st.spinner("モデルを構築中..."):
            # 骨格の取得
            from rdkit.Chem.Scaffolds import MurckoScaffold
            from collections import defaultdict
            
            def get_scaffold(smiles):
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return None
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                return Chem.MolToSmiles(scaffold)
            
            # 骨格を計算
            df['scaffold'] = df['SMILES'].apply(get_scaffold)
            df_clean = df.dropna(subset=descriptor_names + ['ΔΔG.expt.', 'scaffold']).copy()
            
            # 骨格ごとに分類（ユニーク骨格 = 1件のみ）
            scaffold_groups = defaultdict(list)
            for idx, scaffold in df_clean['scaffold'].items():
                scaffold_groups[scaffold].append(idx)
            
            singleton_scaffolds = [s for s, ids in scaffold_groups.items() if len(ids) == 1]
            unique_indices = [scaffold_groups[s][0] for s in singleton_scaffolds]
            non_unique_indices = df_clean.index.difference(unique_indices)
            
            df_train = df_clean.loc[non_unique_indices]
            df_test = df_clean.loc[unique_indices]
            
            X_train = df_train[descriptor_names]
            y_train = df_train['ΔΔG.expt.']
            X_test = df_test[descriptor_names]
            y_test = df_test['ΔΔG.expt.']
            
            # モデル構築
            rf_model = RandomForestRegressor(random_state=42)
            rf_model.fit(X_train, y_train)
            y_pred_train = rf_model.predict(X_train)
            y_pred_test = rf_model.predict(X_test)
        
        # 評価指標の表示
        st.subheader("評価指標")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("テストデータ MSE", f"{mean_squared_error(y_test, y_pred_test):.3f}")
            st.metric("テストデータ R²", f"{r2_score(y_test, y_pred_test):.3f}")
        
        with col2:
            st.metric("トレーニングデータ MSE", f"{mean_squared_error(y_train, y_pred_train):.3f}")
            st.metric("トレーニングデータ R²", f"{r2_score(y_train, y_pred_train):.3f}")
        
        # 分割データの統計情報
        st.subheader("分割データの統計情報")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**テストデータ (ユニーク骨格)**")
            st.metric("サンプル数", len(df_test))
            st.metric("ΔΔG.expt. 範囲", f"{y_test.min():.3f} ~ {y_test.max():.3f}")
            st.metric("ΔΔG.expt. 平均", f"{y_test.mean():.3f}")
            st.metric("ユニーク骨格数", len(set(df_test['scaffold'])))
        
        with col2:
            st.write("**学習データ (複数例がある骨格)**")
            st.metric("サンプル数", len(df_train))
            st.metric("ΔΔG.expt. 範囲", f"{y_train.min():.3f} ~ {y_train.max():.3f}")
            st.metric("ΔΔG.expt. 平均", f"{y_train.mean():.3f}")
            st.metric("ユニーク骨格数", len(set(df_train['scaffold'])))
        
        # モデル情報
        st.subheader("モデル情報")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**使用したモデル**")
            st.write("Random Forest Regressor")
            st.write("**パラメータ**")
            st.write(f"- n_estimators: {rf_model.n_estimators}")
            st.write(f"- max_depth: {rf_model.max_depth}")
            st.write(f"- min_samples_split: {rf_model.min_samples_split}")
            st.write(f"- min_samples_leaf: {rf_model.min_samples_leaf}")
            st.write(f"- random_state: {rf_model.random_state}")
        
        with col2:
            st.write("**データ分割情報**")
            st.write("骨格による分割")
            st.write(f"- 学習データ: {len(X_train)} サンプル (複数例がある骨格)")
            st.write(f"- テストデータ: {len(X_test)} サンプル (ユニーク骨格)")
            st.write(f"- 使用記述子数: {len(descriptor_names)}")
            st.write(f"- 分割基準: Murcko骨格")
        
        # 可視化
        st.subheader("可視化")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # 実測値 vs 予測値
        ax1.scatter(y_train, y_pred_train, label='Training', alpha=0.6)
        ax1.scatter(y_test, y_pred_test, label='Unique Scaffold Test', alpha=0.8)
        min_val = min(y_train.min(), y_test.min())
        max_val = max(y_train.max(), y_test.max())
        ax1.plot([min_val, max_val], [min_val, max_val], 'k--')
        ax1.set_xlabel('Actual ΔΔG.expt.')
        ax1.set_ylabel('Predicted ΔΔG.expt.')
        ax1.set_title('Actual vs Predicted (Scaffold Split)')
        ax1.legend()
        ax1.grid(True)
        
        # 残差プロット
        residuals = y_test - y_pred_test
        ax2.scatter(y_pred_test, residuals)
        ax2.axhline(0, color='black', linestyle='--')
        ax2.set_xlabel('Predicted ΔΔG.expt.')
        ax2.set_ylabel('Residuals')
        ax2.set_title('Residual Plot (Unique Scaffold)')
        ax2.grid(True)
        
        plt.tight_layout()
        st.pyplot(fig)
        
        # 特徴量重要度
        st.subheader("特徴量重要度")
        feature_importance = pd.DataFrame({
            'Feature': descriptor_names,
            'Importance': rf_model.feature_importances_
        }).sort_values('Importance', ascending=False)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(feature_importance['Feature'], feature_importance['Importance'])
        ax.set_xlabel('Importance')
        ax.set_title('Feature Importance')
        plt.tight_layout()
        st.pyplot(fig)
        
        # 詳細な結果
        st.subheader("詳細な結果")
        results_df = pd.DataFrame({
            'Index': X_test.index,
            'SMILES': df_test['SMILES'],
            'Scaffold': df_test['scaffold'],
            'Actual': y_test,
            'Predicted': y_pred_test,
            'Residual': y_test - y_pred_test,
            'Abs_Residual': np.abs(y_test - y_pred_test)
        })
        st.dataframe(results_df)
        
        # 外挿性能の評価
        st.subheader("外挿性能の評価")
        st.write("新しい骨格に対する予測性能:")
        
        # 予測誤差の分析
        mae = np.mean(np.abs(y_test - y_pred_test))
        st.metric("平均絶対誤差 (MAE)", f"{mae:.3f}")
        
        # 予測精度の分布
        accuracy_within_1 = np.mean(np.abs(y_test - y_pred_test) < 1.0) * 100
        accuracy_within_05 = np.mean(np.abs(y_test - y_pred_test) < 0.5) * 100
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric("誤差±1.0以内の予測精度", f"{accuracy_within_1:.1f}%")
        with col2:
            st.metric("誤差±0.5以内の予測精度", f"{accuracy_within_05:.1f}%")
        
        # 骨格の分析
        st.subheader("骨格の分析")
        st.write("**最も予測困難な骨格 (誤差の大きい順)**")
        error_analysis = results_df.nlargest(5, 'Abs_Residual')
        
        for i, (idx, row) in enumerate(error_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                row['Scaffold'], 
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()
        
        st.write("**最も予測精度の高い骨格 (誤差の小さい順)**")
        accuracy_analysis = results_df.nsmallest(5, 'Abs_Residual')
        
        for i, (idx, row) in enumerate(accuracy_analysis.iterrows()):
            st.write(f"**第{i+1}位 (残差: {row['Residual']:.3f})**")
            display_molecule_with_info(
                row['SMILES'], 
                row['Scaffold'], 
                row['Actual'], 
                row['Predicted'], 
                row['Residual']
            )
            st.divider()

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "InterpolationExtrapolation"  # 正しいカテゴリーキーを指定
        
    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()