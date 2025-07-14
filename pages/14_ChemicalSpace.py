import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64
from io import BytesIO

from logic.mol_loader import MoleculeDataLoader

# 共通のPCAプロット作成関数
def create_pca_plot(sample_df, pca_data, pca_model, title_prefix, feature_name):
    """PCAプロットを作成する共通関数"""
    # PCA座標をデータフレームに追加
    sample_df[f'{feature_name}_PC1'] = pca_data[:, 0]
    sample_df[f'{feature_name}_PC2'] = pca_data[:, 1]
    
    fig = go.Figure()
    
    sources = sample_df['Source'].unique()
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, source in enumerate(sources):
        source_data = sample_df[sample_df['Source'] == source]
        
        hover_template = (
            "<b>%{customdata[0]}</b><br>"
            "SMILES: %{customdata[1]}<br>"
            f"{feature_name} PC1: %{{x:.2f}}<br>"
            f"{feature_name} PC2: %{{y:.2f}}<br>"
            "分子量: %{customdata[2]:.2f} Da<br>"
            "LogP: %{customdata[3]:.2f}<br>"
            "TPSA: %{customdata[4]:.2f} Ų<br>"
            "<extra></extra>"
        )
        
        # 分子量に基づくサイズを計算（6-15の範囲にクリップ）
        marker_sizes = np.clip(source_data['MW']/25, 6, 15)
        
        fig.add_trace(go.Scatter(
            x=source_data[f'{feature_name}_PC1'],
            y=source_data[f'{feature_name}_PC2'],
            mode='markers',
            name=source,
            marker=dict(
                color=colors[i % len(colors)],
                size=marker_sizes,
                opacity=0.7,
                sizemin=6,
                line=dict(width=1, color='white')
            ),
            customdata=source_data[['Source', 'SMILES', 'MW', 'LogP', 'TPSA']].values,
            hovertemplate=hover_template
        ))
    
    fig.update_layout(
        title=f'{title_prefix}（PCA）- 寄与率: PC1={pca_model.explained_variance_ratio_[0]:.1%}, PC2={pca_model.explained_variance_ratio_[1]:.1%}',
        xaxis_title=f'{feature_name} PC1 ({pca_model.explained_variance_ratio_[0]:.1%})',
        yaxis_title=f'{feature_name} PC2 ({pca_model.explained_variance_ratio_[1]:.1%})',
        width=800,
        height=600,
        hovermode='closest',
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )
    
    return fig

def display_pca_statistics(pca_model, valid_sample_df, feature_name):
    """PCA統計情報を表示する共通関数"""
    st.subheader(f"📈 {feature_name}記述子のPCA統計情報")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("処理した分子数", len(valid_sample_df))
        st.metric(f"{feature_name}次元数", pca_model.n_components_)
    
    with col2:
        st.metric("PC1寄与率", f"{pca_model.explained_variance_ratio_[0]:.1%}")
        st.metric("PC2寄与率", f"{pca_model.explained_variance_ratio_[1]:.1%}")
    
    with col3:
        st.metric("累積寄与率", f"{pca_model.explained_variance_ratio_[:2].sum():.1%}")
        st.metric("データソース数", len(valid_sample_df['Source'].unique()))

def display_molecule_details(sample_df, feature_df, feature_name, key_prefix):
    """分子詳細情報を表示する共通関数"""
    st.subheader("🔍 選択された分子の詳細情報")
    
    # データソースフィルタリング
    source_filter = st.selectbox(
        "データソースでフィルタ:",
        ["全て"] + list(sample_df['Source'].unique()),
        key=f"{key_prefix}_source_filter"
    )
    
    # フィルタリングを適用
    filtered_df = sample_df.copy()
    
    if source_filter != "全て":
        filtered_df = filtered_df[filtered_df['Source'] == source_filter]
    
    st.write(f"フィルタ後の分子数: {len(filtered_df)} / {len(sample_df)}")
    
    if len(filtered_df) > 0:
        # 分子を選択するためのセレクトボックス
        available_ids = filtered_df.index.tolist()
        
        selected_id = st.selectbox(
            f"分子IDを選択して構造を確認 ({len(available_ids)} 件表示):", 
            available_ids,
            format_func=lambda x: f"ID: {x} - {filtered_df.loc[x, 'Source']} - MW: {filtered_df.loc[x, 'MW']:.1f} - SMILES: {filtered_df.loc[x, 'SMILES'][:30]}...",
            key=f"{key_prefix}_molecule_selector"
        )
        
        # 選択された分子の情報を表示
        selected_mol = filtered_df.loc[selected_id]
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.write("**基本情報:**")
            st.write(f"- データソース: {selected_mol['Source']}")
            st.write(f"- SMILES: `{selected_mol['SMILES']}`")
            st.write(f"- 分子量: {selected_mol['MW']:.2f} Da")
            st.write(f"- LogP: {selected_mol['LogP']:.2f}")
            st.write(f"- TPSA: {selected_mol['TPSA']:.2f} Ų")
            st.write(f"- HBD: {selected_mol['HBD']}")
            st.write(f"- HBA: {selected_mol['HBA']}")
            
            st.write("**PCA座標:**")
            st.write(f"- {feature_name} PC1: {selected_mol[f'{feature_name}_PC1']:.3f}")
            st.write(f"- {feature_name} PC2: {selected_mol[f'{feature_name}_PC2']:.3f}")
        
        with col2:
            st.write("**分子構造:**")
            if selected_mol['Mol'] is not None:
                # RDKitで直接画像を生成してStreamlitで表示
                mol_img = Draw.MolToImage(selected_mol['Mol'], size=(300, 300))
                st.image(mol_img, caption="分子構造", use_container_width=True)
            else:
                st.error("分子構造を生成できませんでした")
        
        return selected_id
    else:
        st.warning("フィルタ条件に一致する分子が見つかりません。")
        return None

# アプリの定義

def mol_to_base64_image(mol, size=(200, 200)):
    """分子をbase64エンコードされた画像に変換"""
    if mol is None:
        return None
    
    try:
        # RDKitで分子画像を生成（PNG形式）
        img_data = Draw.MolToImage(mol, size=size, format='PNG')
        
        # BytesIOを使ってバイト形式に変換
        buffer = BytesIO()
        img_data.save(buffer, format='PNG')
        buffer.seek(0)
        
        # base64エンコード
        img_base64 = base64.b64encode(buffer.getvalue()).decode()
        return f"data:image/png;base64,{img_base64}"
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return None


def load_data():
    file_paths = ['data/TasteDB.smi', 'data/FragranceDB.smi']
    data_frames = []
    
    for i, file in enumerate(file_paths):
        df = MoleculeDataLoader.load_smiles_file(file, has_name=False)
        # データソースを区別するためのラベルを追加
        source_name = "TasteDB" if i == 0 else "FragranceDB"
        df['Source'] = source_name
        
        # 分子記述子を計算
        df['MW'] = df['Mol'].apply(lambda mol: Descriptors.MolWt(mol) if mol else None)
        df['LogP'] = df['Mol'].apply(lambda mol: Descriptors.MolLogP(mol) if mol else None)
        df['TPSA'] = df['Mol'].apply(lambda mol: Descriptors.TPSA(mol) if mol else None)
        df['HBD'] = df['Mol'].apply(lambda mol: Descriptors.NumHDonors(mol) if mol else None)
        df['HBA'] = df['Mol'].apply(lambda mol: Descriptors.NumHAcceptors(mol) if mol else None)
        
        data_frames.append(df)
    
    data = pd.concat(data_frames, ignore_index=True)
    # NaNを含む行を除去
    data = data.dropna(subset=['MW', 'LogP', 'TPSA', 'HBD', 'HBA'])
    return data

def get_common_sample_data(df, sample_size=200, random_state=42):
    """全ての関数で共通のサンプルデータを取得する関数"""
    # 同じランダムシードで同じサンプルを取得
    sample_df = df.sample(n=min(sample_size, len(df)), random_state=random_state).copy()
    return sample_df


def chemical_space_comparison_display():
    st.title("📊 異なる化学空間の比較 (RDKit記述子)")
    st.write("ここでは、異なる化学空間をRDKit記述子で比較し、図示します。")

    df = load_data()

    st.dataframe(df[['SMILES', 'Source', 'MW', 'LogP', 'TPSA', 'HBD', 'HBA']].head(10))

    if len(df) > 0:
        # 共通のサンプルデータを取得
        sample_df = get_common_sample_data(df)
        
        st.write(f"**使用データ:** 全体 {len(df)} 分子から {len(sample_df)} 分子をサンプリング")
        
        # 分子量とLogPの散布図を作成
        st.subheader("📈 RDKit記述子による化学空間分析（分子量 vs LogP）")
        
        fig = go.Figure()
        
        sources = sample_df['Source'].unique()
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        
        for i, source in enumerate(sources):
            source_data = sample_df[sample_df['Source'] == source]
            
            hover_template = (
                "<b>%{customdata[0]}</b><br>"
                "SMILES: %{customdata[1]}<br>"
                "分子量: %{x:.2f} Da<br>"
                "LogP: %{y:.2f}<br>"
                "TPSA: %{customdata[2]:.2f} Ų<br>"
                "HBD: %{customdata[3]}<br>"
                "HBA: %{customdata[4]}<br>"
                "<extra></extra>"
            )
            
            # TPSAに基づくサイズを計算（6-15の範囲にクリップ）
            marker_sizes = np.clip(source_data['TPSA']/10, 6, 15)
            
            fig.add_trace(go.Scatter(
                x=source_data['MW'],
                y=source_data['LogP'],
                mode='markers',
                name=source,
                marker=dict(
                    color=colors[i % len(colors)],
                    size=marker_sizes,
                    opacity=0.7,
                    sizemin=6,
                    line=dict(width=1, color='white')
                ),
                customdata=source_data[['Source', 'SMILES', 'TPSA', 'HBD', 'HBA']].values,
                hovertemplate=hover_template
            ))
        
        fig.update_layout(
            title='RDKit記述子化学空間比較（分子量 vs LogP）',
            xaxis_title='分子量 (Da)',
            yaxis_title='LogP',
            width=800,
            height=600,
            hovermode='closest',
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=0.01
            )
        )
        
        # Streamlitでプロット表示
        st.plotly_chart(fig, use_container_width=True)
        
        # 基本統計情報を表示
        st.subheader("📈 RDKit記述子の基本統計情報")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("処理した分子数", len(sample_df))
            st.metric("平均分子量", f"{sample_df['MW'].mean():.1f} Da")
        
        with col2:
            st.metric("平均LogP", f"{sample_df['LogP'].mean():.2f}")
            st.metric("平均TPSA", f"{sample_df['TPSA'].mean():.1f} Ų")
        
        with col3:
            st.metric("データソース数", len(sample_df['Source'].unique()))
            st.metric("平均HBA", f"{sample_df['HBA'].mean():.1f}")
        
        # 分子詳細情報を表示（PCA座標なしバージョン）
        st.subheader("🔍 選択された分子の詳細情報")
        
        # データソースフィルタリング
        source_filter = st.selectbox(
            "データソースでフィルタ:",
            ["全て"] + list(sample_df['Source'].unique()),
            key="rdkit_source_filter"
        )
        
        # フィルタリングを適用
        filtered_df = sample_df.copy()
        
        if source_filter != "全て":
            filtered_df = filtered_df[filtered_df['Source'] == source_filter]
        
        st.write(f"フィルタ後の分子数: {len(filtered_df)} / {len(sample_df)}")
        
        if len(filtered_df) > 0:
            # 分子を選択するためのセレクトボックス
            available_ids = filtered_df.index.tolist()
            
            selected_id = st.selectbox(
                f"分子IDを選択して構造を確認 ({len(available_ids)} 件表示):", 
                available_ids,
                format_func=lambda x: f"ID: {x} - {filtered_df.loc[x, 'Source']} - MW: {filtered_df.loc[x, 'MW']:.1f} - SMILES: {filtered_df.loc[x, 'SMILES'][:30]}...",
                key="rdkit_molecule_selector"
            )
            
            # 選択された分子の情報を表示
            selected_mol = filtered_df.loc[selected_id]
            
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.write("**基本情報:**")
                st.write(f"- データソース: {selected_mol['Source']}")
                st.write(f"- SMILES: `{selected_mol['SMILES']}`")
                st.write(f"- 分子量: {selected_mol['MW']:.2f} Da")
                st.write(f"- LogP: {selected_mol['LogP']:.2f}")
                st.write(f"- TPSA: {selected_mol['TPSA']:.2f} Ų")
                st.write(f"- HBD: {selected_mol['HBD']}")
                st.write(f"- HBA: {selected_mol['HBA']}")
            
            with col2:
                st.write("**分子構造:**")
                if selected_mol['Mol'] is not None:
                    # RDKitで直接画像を生成してStreamlitで表示
                    mol_img = Draw.MolToImage(selected_mol['Mol'], size=(300, 300))
                    st.image(mol_img, caption="分子構造", use_container_width=True)
                else:
                    st.error("分子構造を生成できませんでした")
        
        else:
            st.warning("フィルタ条件に一致する分子が見つかりません。")
        
        # 複数分子の構造を一度に表示
        st.subheader("🧪 ランダム分子ギャラリー")
        
        # セッション状態でランダムシードを管理
        if 'random_seed' not in st.session_state:
            st.session_state.random_seed = 123
        
        # 新しい分子を表示するボタン
        if st.button("🔄 新しい分子を表示", key="refresh_gallery"):
            st.session_state.random_seed += 1
        
        # ランダムに5つの分子を選択（現在のシードを使用）
        random_mols = sample_df.sample(n=min(5, len(sample_df)), random_state=st.session_state.random_seed)
        
        cols = st.columns(5)
        for i, (mol_id, mol_data) in enumerate(random_mols.iterrows()):
            with cols[i]:
                if mol_data['Mol'] is not None:
                    mol_img = Draw.MolToImage(mol_data['Mol'], size=(150, 150))
                    st.image(mol_img, caption=f"ID: {mol_id}\n{mol_data['Source']}\nMW: {mol_data['MW']:.1f}")
                    st.caption(f"LogP: {mol_data['LogP']:.2f}")
                else:
                    st.error(f"ID: {mol_id}\n構造生成エラー")
        
        st.info("💡 プロット上の点にマウスを重ねると詳細情報が表示されます。点のサイズはTPSAを表しています。")
        
    else:
        st.warning("データが見つかりません。")


def chemical_space_mqn_comparison_display():
    st.title("📊 異なる化学空間の比較 (MQN記述子)")
    st.write("ここでは、異なる化学空間をMQN（Molecular Quantum Numbers）で比較し、図示します。")

    try:
        from rdkit.Chem import rdMolDescriptors
        mqn_available = True
    except ImportError:
        mqn_available = False
        st.error("MQNモジュールがインストールされていません。")

    df = load_data()
    
    if len(df) > 0 and mqn_available:
        # 共通のサンプルデータを取得
        sample_df = get_common_sample_data(df)
        
        st.write(f"**使用データ:** 全体 {len(df)} 分子から {len(sample_df)} 分子をサンプリング")
        
        # MQN記述子を計算
        st.info("MQN記述子を計算中...")
        mqn_vectors = []
        valid_indices = []
        
        for idx, row in sample_df.iterrows():
            mol = row['Mol']
            if mol is not None:
                try:
                    mqn = rdMolDescriptors.MQNs_(mol)
                    mqn_list = list(mqn)
                    if len(mqn_list) == 42:  # MQNは42次元
                        mqn_vectors.append(mqn_list)
                        valid_indices.append(idx)
                except Exception as e:
                    continue
        
        if len(mqn_vectors) > 0:
            import numpy as np
            from sklearn.decomposition import PCA
            
            # MQNベクトルをDataFrameに変換
            mqn_df = pd.DataFrame(mqn_vectors, index=valid_indices)
            mqn_df.columns = [f'MQN_{i}' for i in range(42)]
            
            # 有効なデータのみを取得
            valid_sample_df = sample_df.loc[valid_indices].copy()
            
            # PCAで2次元に次元削減
            pca = PCA(n_components=2)
            mqn_pca = pca.fit_transform(mqn_vectors)
            
            # 共通関数を使用してPCAプロットを作成
            st.subheader("📈 MQN記述子による化学空間分析（PCA）")
            fig = create_pca_plot(valid_sample_df, mqn_pca, pca, "MQN記述子化学空間比較", "MQN")
            st.plotly_chart(fig, use_container_width=True)
            
            # PCA統計情報を表示
            display_pca_statistics(pca, valid_sample_df, "MQN")
            
            # 主成分負荷量の分析
            st.subheader("🎯 MQN記述子の寄与度分析")
            
            # MQN記述子の名前を定義（42次元）
            mqn_names = [
                'MW_100', 'MW_200', 'MW_300', 'MW_400', 'MW_500', 'MW_600', 'MW_700', 'MW_800', 'MW_900',
                'AMW', 'N_atoms', 'N_bonds', 'N_BM', 'N_BO', 'N_SC', 'N_DB', 'N_AR', 'N_TB', 'N_carb',
                'N_het', 'PC_1', 'PC_2', 'PC_3', 'PC_4', 'PC_5', 'PC_6', 'A_acc', 'A_don', 'A_acid',
                'A_base', 'A_ampho', 'A_aro6', 'A_aro7plus', 'A_sat', 'A_unsat', 'CS', 'HetCarbonRatio',
                'RPCG', 'C1SP1', 'C2SP1', 'C1SP2', 'C2SP2'
            ]
            
            # 主成分負荷量（特徴量の重み）を取得
            pc1_loadings = pca.components_[0]
            pc2_loadings = pca.components_[1]
            
            # 寄与度の絶対値でソート
            pc1_importance = pd.DataFrame({
                'MQN記述子': mqn_names,
                '負荷量': pc1_loadings,
                '絶対値': np.abs(pc1_loadings)
            }).sort_values('絶対値', ascending=False)
            
            pc2_importance = pd.DataFrame({
                'MQN記述子': mqn_names,
                '負荷量': pc2_loadings,
                '絶対値': np.abs(pc2_loadings)
            }).sort_values('絶対値', ascending=False)
            
            # 寄与度上位を表示
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**PC1への寄与度上位10位**")
                top10_pc1 = pc1_importance.head(10)[['MQN記述子', '負荷量']].copy()
                top10_pc1['負荷量'] = top10_pc1['負荷量'].round(3)
                st.dataframe(top10_pc1, use_container_width=True, hide_index=True)
            
            with col2:
                st.write("**PC2への寄与度上位10位**")
                top10_pc2 = pc2_importance.head(10)[['MQN記述子', '負荷量']].copy()
                top10_pc2['負荷量'] = top10_pc2['負荷量'].round(3)
                st.dataframe(top10_pc2, use_container_width=True, hide_index=True)
            
            # 分子詳細情報を表示
            selected_id = display_molecule_details(valid_sample_df, mqn_df, "MQN", "mqn")
            
            # 選択された分子のMQN詳細情報（オプショナル）
            if selected_id is not None:
                selected_mqn_features = mqn_df.loc[selected_id]
                
                # MQN記述子の詳細表示
                with st.expander("📊 選択された分子のMQN記述子詳細"):
                    mqn_col1, mqn_col2, mqn_col3 = st.columns(3)
                    
                    with mqn_col1:
                        st.write("**原子・結合数:**")
                        st.write(f"- 原子数: {selected_mqn_features['MQN_10']}")
                        st.write(f"- 結合数: {selected_mqn_features['MQN_11']}")
                        st.write(f"- 炭素原子数: {selected_mqn_features['MQN_18']}")
                        st.write(f"- ヘテロ原子数: {selected_mqn_features['MQN_19']}")
                    
                    with mqn_col2:
                        st.write("**化学的性質:**")
                        st.write(f"- H結合アクセプター: {selected_mqn_features['MQN_26']}")
                        st.write(f"- H結合ドナー: {selected_mqn_features['MQN_27']}")
                        st.write(f"- 酸性基数: {selected_mqn_features['MQN_28']}")
                        st.write(f"- 塩基性基数: {selected_mqn_features['MQN_29']}")
                    
                    with mqn_col3:
                        st.write("**環構造:**")
                        st.write(f"- 6員芳香環: {selected_mqn_features['MQN_31']}")
                        st.write(f"- 7員環以上: {selected_mqn_features['MQN_32']}")
                        st.write(f"- 飽和環: {selected_mqn_features['MQN_33']}")
                        st.write(f"- 不飽和環: {selected_mqn_features['MQN_34']}")
            
            st.info("💡 MQN記述子は分子の構造特徴を42次元の数値ベクトルで表現します。点のサイズは分子量を表しています。")
            
        else:
            st.error("有効なMQN記述子を計算できませんでした。")
    
    elif not mqn_available:
        st.warning("MQNモジュールが利用できません。")
    else:
        st.warning("データが見つかりません。")


def chemical_space_smifp_comparison_display():
    st.title("📊 異なる化学空間の比較 (SELFIES記述子)")
    st.write("ここでは、SMILESをSELFIESに変換してから特徴量抽出を行い、化学空間を比較・図示します。")

    # SELFIES対応表データ
    data = [
        {"no": 1, "SMILES symbol": "C", "feature counted": "非芳香族C", "SELFIESでのカウント可否": "✅ [C]", "備考": ""},
        {"no": 2, "SMILES symbol": "c", "feature counted": "芳香族C", "SELFIESでのカウント可否": "✅ [c]", "備考": ""},
        {"no": 3, "SMILES symbol": "N", "feature counted": "非芳香族N", "SELFIESでのカウント可否": "✅ [N]", "備考": ""},
        {"no": 4, "SMILES symbol": "n", "feature counted": "芳香族N", "SELFIESでのカウント可否": "✅ [n]", "備考": ""},
        {"no": 5, "SMILES symbol": "O", "feature counted": "非芳香族O", "SELFIESでのカウント可否": "✅ [O]", "備考": ""},
        {"no": 6, "SMILES symbol": "o", "feature counted": "芳香族O", "SELFIESでのカウント可否": "✅ [o]", "備考": ""},
        {"no": 7, "SMILES symbol": "S", "feature counted": "非芳香族S", "SELFIESでのカウント可否": "✅ [S]", "備考": ""},
        {"no": 8, "SMILES symbol": "s", "feature counted": "芳香族S", "SELFIESでのカウント可否": "✅ [s]", "備考": ""},
        {"no": 9, "SMILES symbol": "F", "feature counted": "F原子", "SELFIESでのカウント可否": "✅ [F]", "備考": ""},
        {"no": 10, "SMILES symbol": "Cl", "feature counted": "Cl原子", "SELFIESでのカウント可否": "✅ [Cl]", "備考": ""},
        {"no": 11, "SMILES symbol": "Br", "feature counted": "Br原子", "SELFIESでのカウント可否": "✅ [Br]", "備考": ""},
        {"no": 12, "SMILES symbol": "I", "feature counted": "I原子", "SELFIESでのカウント可否": "✅ [I]", "備考": ""},
        {"no": 13, "SMILES symbol": "P", "feature counted": "非芳香族P", "SELFIESでのカウント可否": "✅ [P]", "備考": ""},
        {"no": 14, "SMILES symbol": "p", "feature counted": "芳香族P", "SELFIESでのカウント可否": "✅ [p]", "備考": ""},
        {"no": 15, "SMILES symbol": "B", "feature counted": "B原子", "SELFIESでのカウント可否": "✅ [B]", "備考": ""},
        {"no": 16, "SMILES symbol": "X", "feature counted": "その他の原子", "SELFIESでのカウント可否": "⭕ 未知トークンとして処理", "備考": ""},
        {"no": 17, "SMILES symbol": "-", "feature counted": "単結合", "SELFIESでのカウント可否": "省略", "備考": "SELFIESでは結合を明示しない"},
        {"no": 18, "SMILES symbol": "=", "feature counted": "二重結合", "SELFIESでのカウント可否": "✅ [=C] など", "備考": "`=` を含むトークンを集計"},
        {"no": 19, "SMILES symbol": "#", "feature counted": "三重結合", "SELFIESでのカウント可否": "✅ [#C] など", "備考": "`#` を含むトークンを集計"},
        {"no": 20, "SMILES symbol": "[", "feature counted": "特殊記号", "SELFIESでのカウント可否": "省略", "備考": "SELFIESでは使わない"},
        {"no": 21, "SMILES symbol": "-", "feature counted": "負電荷", "SELFIESでのカウント可否": "✅ [Charge-1] など", "備考": "`-1` を含むトークン"},
        {"no": 22, "SMILES symbol": "+", "feature counted": "正電荷", "SELFIESでのカウント可否": "✅ [Charge+1] など", "備考": "`+1` を含むトークン"},
        {"no": 23, "SMILES symbol": "H", "feature counted": "明示的水素", "SELFIESでのカウント可否": "✅ [H]", "備考": ""},
        {"no": 24, "SMILES symbol": "(", "feature counted": "分岐点", "SELFIESでのカウント可否": "✅ [Branch1], [Branch2]", "備考": ""},
        {"no": "25–34", "SMILES symbol": "1〜%", "feature counted": "環のサイズ情報", "SELFIESでのカウント可否": "🔶 [Ring1]〜[Ring9]", "備考": "環の存在/数は可、サイズは不可"},
    ]

    # SELFIES対応表の表示
    df_selfies = pd.DataFrame(data)
    
    with st.expander("📋 SMIfp項目とSELFIESトークンの対応表を表示"):
        st.write("以下は、SMIfp項目とSELFIESトークンの対応表です。")
        st.dataframe(df_selfies, use_container_width=True)

    # データの読み込み
    df = load_data()

    # SELFIES特徴量抽出の確認
    try:
        from logic.selfies_features import count_selfies_features
        selfies_available = True
        st.success("SELFIES特徴量抽出ライブラリが利用可能です")
    except ImportError:
        selfies_available = False
        st.error("SELFIES特徴量抽出ライブラリが見つかりません")
        return

    if len(df) > 0 and selfies_available:
        # 共通のサンプルデータを取得
        sample_df = get_common_sample_data(df)
        
        st.write(f"**使用データ:** 全体 {len(df)} 分子から {len(sample_df)} 分子をサンプリング")
        
        # SELFIES特徴量を計算
        st.info("SELFIES特徴量を計算中...")
        selfies_features_list = []
        valid_indices = []
        
        progress_bar = st.progress(0)
        total_molecules = len(sample_df)
        
        for i, (idx, row) in enumerate(sample_df.iterrows()):
            smiles = row['SMILES']
            if smiles:
                try:
                    features = count_selfies_features(smiles)
                    if features is not None:
                        selfies_features_list.append(features)
                        valid_indices.append(idx)
                except Exception as e:
                    continue
            
            # プログレスバーを更新
            progress_bar.progress((i + 1) / total_molecules)
        
        progress_bar.empty()
        
        if len(selfies_features_list) > 0:
            import numpy as np
            from sklearn.decomposition import PCA
            
            # SELFIES特徴量をDataFrameに変換
            selfies_df = pd.DataFrame(selfies_features_list, index=valid_indices)
            
            # 有効なデータのみを取得
            valid_sample_df = sample_df.loc[valid_indices].copy()
            
            # PCAで2次元に次元削減
            features_for_pca = selfies_df.select_dtypes(include=[np.number]).fillna(0)
            
            if features_for_pca.shape[1] >= 2:
                pca = PCA(n_components=2)
                selfies_pca = pca.fit_transform(features_for_pca)
                
                # 共通関数を使用してPCAプロットを作成
                st.subheader("📈 SELFIES特徴量による化学空間分析（PCA）")
                fig = create_pca_plot(valid_sample_df, selfies_pca, pca, "SELFIES特徴量化学空間比較", "SELFIES")
                st.plotly_chart(fig, use_container_width=True)
                
                # PCA統計情報を表示
                display_pca_statistics(pca, valid_sample_df, "SELFIES")
                
                # 分子詳細情報を表示
                selected_id = display_molecule_details(valid_sample_df, selfies_df, "SELFIES", "selfies")
                
                # 選択された分子のSELFIES詳細情報（オプショナル）
                if selected_id is not None:
                    selected_selfies_features = selfies_df.loc[selected_id]
                    
                    # SELFIES特徴量の詳細表示
                    with st.expander("📊 選択された分子のSELFIES特徴量詳細"):
                        selfies_col1, selfies_col2 = st.columns(2)
                        
                        with selfies_col1:
                            st.write("**原子数:**")
                            if 'C_aliphatic' in selected_selfies_features:
                                st.write(f"- 脂肪族C: {selected_selfies_features['C_aliphatic']}")
                            if 'C_aromatic' in selected_selfies_features:
                                st.write(f"- 芳香族C: {selected_selfies_features['C_aromatic']}")
                            if 'N_aliphatic' in selected_selfies_features:
                                st.write(f"- 脂肪族N: {selected_selfies_features['N_aliphatic']}")
                            if 'N_aromatic' in selected_selfies_features:
                                st.write(f"- 芳香族N: {selected_selfies_features['N_aromatic']}")
                            if 'O_aliphatic' in selected_selfies_features:
                                st.write(f"- 脂肪族O: {selected_selfies_features['O_aliphatic']}")
                            if 'O_aromatic' in selected_selfies_features:
                                st.write(f"- 芳香族O: {selected_selfies_features['O_aromatic']}")
                            if 'S' in selected_selfies_features:
                                st.write(f"- S: {selected_selfies_features['S']}")
                            if 'F' in selected_selfies_features:
                                st.write(f"- F: {selected_selfies_features['F']}")
                            if 'Cl' in selected_selfies_features:
                                st.write(f"- Cl: {selected_selfies_features['Cl']}")
                            if 'Br' in selected_selfies_features:
                                st.write(f"- Br: {selected_selfies_features['Br']}")
                        
                        with selfies_col2:
                            st.write("**構造特徴:**")
                            if 'double_bonds' in selected_selfies_features:
                                st.write(f"- 二重結合: {selected_selfies_features['double_bonds']}")
                            if 'triple_bonds' in selected_selfies_features:
                                st.write(f"- 三重結合: {selected_selfies_features['triple_bonds']}")
                            if 'branches' in selected_selfies_features:
                                st.write(f"- 分岐数: {selected_selfies_features['branches']}")
                            if 'rings' in selected_selfies_features:
                                st.write(f"- 環数: {selected_selfies_features['rings']}")
                            if 'stereo_centers' in selected_selfies_features:
                                st.write(f"- 立体中心: {selected_selfies_features['stereo_centers']}")
                            if 'explicit_H' in selected_selfies_features:
                                st.write(f"- 明示的H: {selected_selfies_features['explicit_H']}")
                            if 'positive_charge' in selected_selfies_features:
                                st.write(f"- 正電荷: {selected_selfies_features['positive_charge']}")
                            if 'negative_charge' in selected_selfies_features:
                                st.write(f"- 負電荷: {selected_selfies_features['negative_charge']}")
                            if 'special_atoms' in selected_selfies_features:
                                st.write(f"- 特殊原子: {selected_selfies_features['special_atoms']}")
                
                # SELFIES特徴量の説明
                with st.expander("📖 SELFIES特徴量の詳細説明"):
                    st.write("""
                    **SELFIES特徴量について:**
                    - C_aliphatic/C_aromatic: 脂肪族/芳香族炭素原子数
                    - N_aliphatic/N_aromatic: 脂肪族/芳香族窒素原子数
                    - O_aliphatic/O_aromatic: 脂肪族/芳香族酸素原子数
                    - S, P, F, Cl, Br, I, B: 各種原子数
                    - double_bonds/triple_bonds: 二重結合/三重結合数
                    - positive_charge/negative_charge: 正電荷/負電荷数
                    - explicit_H: 明示的水素数
                    - branches: 分岐数
                    - rings: 環の数
                    - stereo_centers: 立体中心数
                    
                    **SELFIES vs SMIfp:**
                    - SMILESの文字列解析ではなく、SELFIESトークンベースで特徴量を抽出
                    - より安定した分子表現により、ノイズの少ない特徴量を獲得
                    - 分子の構造特徴をより直接的に反映
                    """)
                
                st.info("💡 SELFIES特徴量を使用することで、分子構造の特徴をより詳細に分析できます。点のサイズは分子量を表しています。")
                
            else:
                st.error("PCAに十分な特徴量がありません。")
                
        else:
            st.error("有効なSELFIES特徴量を計算できませんでした。")
    
    elif not selfies_available:
        st.warning("SELFIES特徴量抽出機能が利用できません。")
    else:
        st.warning("データが見つかりません。")


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ChemicalSpace"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()