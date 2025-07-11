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


def chemical_space_comparison_display():
    st.title("📊 異なる化学空間の比較")
    st.write("ここでは、異なる化学空間をRDKit descriptorで比較し、図示します。")

    df = load_data()

    st.dataframe(df[['SMILES', 'Source', 'MW', 'LogP', 'TPSA', 'HBD', 'HBA']].head(10))

    # インタラクティブな散布図を表示
    st.subheader("分子量とLogPの散布図（インタラクティブ）")
    
    if len(df) > 0:
        # 処理を軽くするためサンプルデータに限定
        sample_df = df.sample(n=min(200, len(df)), random_state=42).copy().reset_index(drop=True)
        sample_df = df
        
        # Plotlyを使用してインタラクティブな散布図を作成（画像なしバージョン）
        fig = go.Figure()
        
        # データソース別に色分けして追加
        sources = sample_df['Source'].unique()
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        
        for i, source in enumerate(sources):
            source_data = sample_df[sample_df['Source'] == source]
            
            # hoverテンプレートを作成（テキスト情報のみ）
            hover_template = (
                "<b>%{customdata[0]}</b><br>"
                "SMILES: %{customdata[1]}<br>"
                "分子量: %{x:.2f} Da<br>"
                "LogP: %{y:.2f}<br>"
                "TPSA: %{customdata[2]:.2f} Ų<br>"
                "HBD: %{customdata[3]}<br>"
                "HBA: %{customdata[4]}<br>"
                "インデックス: %{customdata[5]}"
                "<extra></extra>"
            )
            
            fig.add_trace(go.Scatter(
                x=source_data['MW'],
                y=source_data['LogP'],
                mode='markers',
                name=source,
                marker=dict(
                    color=colors[i % len(colors)],
                    size=8,
                    opacity=0.7
                ),
                customdata=source_data[['Source', 'SMILES', 'TPSA', 'HBD', 'HBA']].reset_index().values,
                hovertemplate=hover_template
            ))
        
        # グラフの見た目を調整
        fig.update_layout(
            title='分子量 vs LogP（クリックで詳細表示）',
            xaxis_title='分子量 (Da)',
            yaxis_title='LogP',
            width=800,
            height=600,
            hovermode='closest'
        )
        
        # Streamlitでプロット表示
        selected_points = st.plotly_chart(fig, use_container_width=True, on_select="rerun")
        
        # クリックされた点の情報を表示
        st.subheader("🔍 選択された分子の詳細情報")
        
        # データソースフィルタリング
        source_filter = st.selectbox(
            "データソースでフィルタ:",
            ["全て"] + list(sample_df['Source'].unique()),
            key="chemical_space_source_filter"
        )
        
        # フィルタリングを適用
        filtered_df = sample_df.copy()
        
        if source_filter != "全て":
            filtered_df = filtered_df[filtered_df['Source'] == source_filter]
        
        st.write(f"フィルタ後の分子数: {len(filtered_df)} / {len(sample_df)}")
        
        if len(filtered_df) > 0:
            # 分子を選択するためのセレクトボックス（全フィルタ結果から選択可能）
            available_ids = filtered_df.index.tolist()
            
            selected_id = st.selectbox(
                f"分子IDを選択して構造を確認 ({len(available_ids)} 件表示):", 
                available_ids,
                format_func=lambda x: f"ID: {x} - {filtered_df.loc[x, 'Source']} - MW: {filtered_df.loc[x, 'MW']:.1f} - SMILES: {filtered_df.loc[x, 'SMILES'][:30]}...",
                key="chemical_space_molecule_selector"
            )
        else:
            st.warning("フィルタ条件に一致する分子が見つかりません。")
            selected_id = sample_df.index[0]  # デフォルトとして最初の分子を選択
        
        # 選択された分子の情報を表示
        selected_mol = sample_df.loc[selected_id]
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.write("**分子情報:**")
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
        
        st.info("💡 プロット上の点にマウスを重ねると詳細情報が表示されます。分子構造は上の選択ボックスやギャラリーで確認できます。")
        
    else:
        st.warning("データが見つかりません。")


def chemical_space_mqn_comparison_display():
    st.title("📊 異なる化学空間の比較 (MQN)")
    st.write("ここでは、異なる化学空間をMQN（Molecular Quantum Numbers）で比較し、図示します。")

    try:
        from rdkit.Chem import rdMolDescriptors
        mqn_available = True
    except ImportError:
        mqn_available = False
        st.error("MQNモジュールがインストールされていません。")

    df = load_data()
    
    if len(df) > 0 and mqn_available:
        # 処理を軽くするためサンプルデータに限定
        sample_df = df.sample(n=min(200, len(df)), random_state=42).copy()
        
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
            valid_sample_df['mol_image'] = valid_sample_df['Mol'].apply(mol_to_base64_image)
            
            # PCAで2次元に次元削減
            pca = PCA(n_components=2)
            mqn_pca = pca.fit_transform(mqn_vectors)
            
            valid_sample_df['MQN_PC1'] = mqn_pca[:, 0]
            valid_sample_df['MQN_PC2'] = mqn_pca[:, 1]
            
            # Plotlyでの可視化
            fig = go.Figure()
            
            sources = valid_sample_df['Source'].unique()
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
            
            for i, source in enumerate(sources):
                source_data = valid_sample_df[valid_sample_df['Source'] == source]
                
                hover_template = (
                    "<b>%{customdata[0]}</b><br>"
                    "SMILES: %{customdata[1]}<br>"
                    "MQN PC1: %{x:.2f}<br>"
                    "MQN PC2: %{y:.2f}<br>"
                    "分子量: %{customdata[2]:.2f} Da<br>"
                    "LogP: %{customdata[3]:.2f}<br>"
                    "TPSA: %{customdata[4]:.2f} Ų<br>"
                    "<img src='%{customdata[5]}' width='200'>"
                    "<extra></extra>"
                )
                
                # 分子量に基づくサイズを計算（4-20の範囲にクリップ）
                marker_sizes = np.clip(source_data['MW']/10, 4, 20)
                
                fig.add_trace(go.Scatter(
                    x=source_data['MQN_PC1'],
                    y=source_data['MQN_PC2'],
                    mode='markers',
                    name=source,
                    marker=dict(
                        color=colors[i % len(colors)],
                        size=marker_sizes,
                        opacity=0.7,
                        sizemin=4
                    ),
                    customdata=source_data[['Source', 'SMILES', 'MW', 'LogP', 'TPSA', 'mol_image']].values,
                    hovertemplate=hover_template
                ))
            
            fig.update_layout(
                title=f'MQN化学空間比較（PCA）- 寄与率: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}',
                xaxis_title=f'MQN PC1 ({pca.explained_variance_ratio_[0]:.1%})',
                yaxis_title=f'MQN PC2 ({pca.explained_variance_ratio_[1]:.1%})',
                width=800,
                height=600,
                hovermode='closest'
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # MQN統計情報を表示
            st.subheader("📈 MQN記述子の統計情報")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("処理した分子数", len(valid_sample_df))
                st.metric("MQN次元数", 42)
            
            with col2:
                st.metric("PC1寄与率", f"{pca.explained_variance_ratio_[0]:.1%}")
                st.metric("PC2寄与率", f"{pca.explained_variance_ratio_[1]:.1%}")
            
            with col3:
                st.metric("累積寄与率", f"{pca.explained_variance_ratio_[:2].sum():.1%}")
                st.metric("データソース数", len(sources))
            
            # 主成分負荷量（特徴量の重み）の分析
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
                
                # PC1寄与度のバーチャート
                fig_pc1 = go.Figure()
                fig_pc1.add_trace(go.Bar(
                    x=top10_pc1['負荷量'],
                    y=top10_pc1['MQN記述子'],
                    orientation='h',
                    marker_color=['red' if x < 0 else 'blue' for x in top10_pc1['負荷量']],
                    text=[f'{x:.3f}' for x in top10_pc1['負荷量']],
                    textposition='auto'
                ))
                fig_pc1.update_layout(
                    title='PC1への寄与度（上位10位）',
                    xaxis_title='負荷量',
                    yaxis_title='MQN記述子',
                    height=400,
                    yaxis={'categoryorder': 'total ascending'}
                )
                st.plotly_chart(fig_pc1, use_container_width=True)
            
            with col2:
                st.write("**PC2への寄与度上位10位**")
                top10_pc2 = pc2_importance.head(10)[['MQN記述子', '負荷量']].copy()
                top10_pc2['負荷量'] = top10_pc2['負荷量'].round(3)
                st.dataframe(top10_pc2, use_container_width=True, hide_index=True)
                
                # PC2寄与度のバーチャート
                fig_pc2 = go.Figure()
                fig_pc2.add_trace(go.Bar(
                    x=top10_pc2['負荷量'],
                    y=top10_pc2['MQN記述子'],
                    orientation='h',
                    marker_color=['red' if x < 0 else 'blue' for x in top10_pc2['負荷量']],
                    text=[f'{x:.3f}' for x in top10_pc2['負荷量']],
                    textposition='auto'
                ))
                fig_pc2.update_layout(
                    title='PC2への寄与度（上位10位）',
                    xaxis_title='負荷量',
                    yaxis_title='MQN記述子',
                    height=400,
                    yaxis={'categoryorder': 'total ascending'}
                )
                st.plotly_chart(fig_pc2, use_container_width=True)
            
            # MQN記述子の説明
            with st.expander("📖 MQN記述子の詳細説明"):
                st.write("""
                **MQN記述子の各項目説明:**
                
                **分子量関連 (MW_xxx):**
                - MW_100〜MW_900: 分子量が100〜900の範囲にあるかのフラグ
                
                **原子・結合数:**
                - AMW: 平均分子量
                - N_atoms: 原子数
                - N_bonds: 結合数
                - N_BM: 分岐点数
                - N_BO: 結合次数
                - N_SC: 単結合数
                - N_DB: 二重結合数
                - N_AR: 芳香環結合数
                - N_TB: 三重結合数
                - N_carb: 炭素原子数
                - N_het: ヘテロ原子数
                
                **極性・化学的性質:**
                - PC_1〜PC_6: 極性分類
                - A_acc: 水素結合アクセプター数
                - A_don: 水素結合ドナー数
                - A_acid: 酸性基数
                - A_base: 塩基性基数
                - A_ampho: 両性基数
                
                **環構造:**
                - A_aro6: 6員環芳香環数
                - A_aro7plus: 7員環以上の芳香環数
                - A_sat: 飽和環数
                - A_unsat: 不飽和環数
                
                **その他:**
                - CS: カイラル中心数
                - HetCarbonRatio: ヘテロ原子/炭素原子比
                - RPCG: 環形成可能性
                - C1SP1, C2SP1: sp軌道炭素数
                - C1SP2, C2SP2: sp²軌道炭素数
                """)
            
            # 個別のMQN値も表示
            with st.expander("📊 MQN記述子の詳細"):
                st.write("最初の5分子のMQN値:")
                display_mqn = mqn_df.head().round(2)
                st.dataframe(display_mqn)
                
                # MQN記述子の説明
                st.write("""
                **MQN記述子について:**
                - MQN (Molecular Quantum Numbers) は42次元の整数ベクトル
                - 分子の構造的特徴を数値化
                - 原子数、結合数、環の情報などを含む
                - PCAで2次元に圧縮して可視化
                """)
            
            # 選択された分子の詳細情報
            st.subheader("🔍 選択された分子の詳細情報")
            
            # データソースフィルタリング
            source_filter = st.selectbox(
                "データソースでフィルタ:",
                ["全て"] + list(valid_sample_df['Source'].unique()),
                key="mqn_source_filter"
            )
            
            # フィルタリングを適用
            filtered_df = valid_sample_df.copy()
            
            if source_filter != "全て":
                filtered_df = filtered_df[filtered_df['Source'] == source_filter]
            
            st.write(f"フィルタ後の分子数: {len(filtered_df)} / {len(valid_sample_df)}")
            
            if len(filtered_df) > 0:
                # 分子を選択するためのセレクトボックス
                available_ids = filtered_df.index.tolist()
                
                selected_id = st.selectbox(
                    f"分子IDを選択して構造を確認 ({len(available_ids)} 件表示):", 
                    available_ids,
                    format_func=lambda x: f"ID: {x} - {filtered_df.loc[x, 'Source']} - MW: {filtered_df.loc[x, 'MW']:.1f} - SMILES: {filtered_df.loc[x, 'SMILES'][:30]}...",
                    key="mqn_molecule_selector"
                )
                
                # 選択された分子の情報を表示
                selected_mol = filtered_df.loc[selected_id]
                selected_mqn_features = mqn_df.loc[selected_id]
                
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
                    st.write(f"- MQN PC1: {selected_mol['MQN_PC1']:.3f}")
                    st.write(f"- MQN PC2: {selected_mol['MQN_PC2']:.3f}")
                
                with col2:
                    st.write("**分子構造:**")
                    if selected_mol['Mol'] is not None:
                        # RDKitで直接画像を生成してStreamlitで表示
                        mol_img = Draw.MolToImage(selected_mol['Mol'], size=(300, 300))
                        st.image(mol_img, caption="分子構造", use_container_width=True)
                    else:
                        st.error("分子構造を生成できませんでした")
                
                # MQN記述子の詳細表示
                st.write("**MQN記述子詳細:**")
                mqn_col1, mqn_col2, mqn_col3 = st.columns(3)
                
                with mqn_col1:
                    st.write("分子量範囲:")
                    for i in range(9):
                        mqn_key = f'MQN_{i}'
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            if value > 0:
                                st.write(f"- MW_{(i+1)*100}: {value}")
                    
                    st.write("原子・結合数:")
                    mqn_mapping = {
                        'MQN_9': 'AMW (平均分子量)',
                        'MQN_10': 'N_atoms (原子数)',
                        'MQN_11': 'N_bonds (結合数)',
                        'MQN_18': 'N_carb (炭素原子数)',
                        'MQN_19': 'N_het (ヘテロ原子数)'
                    }
                    for mqn_key, description in mqn_mapping.items():
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            st.write(f"- {description}: {value}")
                
                with mqn_col2:
                    st.write("結合タイプ:")
                    bond_mapping = {
                        'MQN_12': 'N_BM (分岐点数)',
                        'MQN_13': 'N_BO (結合次数)',
                        'MQN_14': 'N_SC (単結合数)',
                        'MQN_15': 'N_DB (二重結合数)',
                        'MQN_16': 'N_AR (芳香環結合数)',
                        'MQN_17': 'N_TB (三重結合数)'
                    }
                    for mqn_key, description in bond_mapping.items():
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            if value > 0:
                                st.write(f"- {description}: {value}")
                    
                    st.write("化学的性質:")
                    chem_mapping = {
                        'MQN_26': 'A_acc (H結合アクセプター)',
                        'MQN_27': 'A_don (H結合ドナー)',
                        'MQN_28': 'A_acid (酸性基数)',
                        'MQN_29': 'A_base (塩基性基数)',
                        'MQN_30': 'A_ampho (両性基数)'
                    }
                    for mqn_key, description in chem_mapping.items():
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            if value > 0:
                                st.write(f"- {description}: {value}")
                
                with mqn_col3:
                    st.write("環構造:")
                    ring_mapping = {
                        'MQN_31': 'A_aro6 (6員芳香環)',
                        'MQN_32': 'A_aro7plus (7員環以上)',
                        'MQN_33': 'A_sat (飽和環)',
                        'MQN_34': 'A_unsat (不飽和環)'
                    }
                    for mqn_key, description in ring_mapping.items():
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            if value > 0:
                                st.write(f"- {description}: {value}")
                    
                    st.write("その他:")
                    other_mapping = {
                        'MQN_35': 'CS (カイラル中心数)',
                        'MQN_36': 'HetCarbonRatio',
                        'MQN_37': 'RPCG (環形成可能性)'
                    }
                    for mqn_key, description in other_mapping.items():
                        if mqn_key in selected_mqn_features:
                            value = selected_mqn_features[mqn_key]
                            if value > 0:
                                st.write(f"- {description}: {value:.2f}")
                
                # 全MQN値の表示（オプション）
                with st.expander("📋 全MQN値を表示"):
                    mqn_display_df = pd.DataFrame({
                        'MQN記述子': [f'MQN_{i}' for i in range(42)],
                        '値': [selected_mqn_features[f'MQN_{i}'] for i in range(42)]
                    })
                    st.dataframe(mqn_display_df, use_container_width=True, hide_index=True)
            
            else:
                st.warning("フィルタ条件に一致する分子が見つかりません。")
            
            st.info("💡 マウスオーバーで化学構造とMQN情報を確認できます。点のサイズは分子量を表しています。")
            
        else:
            st.error("有効なMQN記述子を計算できませんでした。")
    
    elif not mqn_available:
        st.warning("MQNモジュールが利用できません。代替として基本記述子を使用します。")
        # フォールバック: 基本記述子を使用
        if len(df) > 0:
            sample_df = df.sample(n=min(150, len(df)), random_state=42).copy()
            fig = px.scatter(
                sample_df,
                x='MW',
                y='HBA',
                color='Source',
                size='TPSA',
                title='MQN風化学空間比較（分子量 vs 水素結合アクセプター数）',
                labels={'MW': '分子量 (Da)', 'HBA': '水素結合アクセプター数'},
                hover_data=['SMILES', 'LogP', 'TPSA', 'HBD']
            )
            fig.update_layout(width=800, height=600, hovermode='closest')
            st.plotly_chart(fig, use_container_width=True)
            st.info("💡 実際のMQN記述子を使用する場合は、適切なライブラリのインストールが必要です。")
    else:
        st.warning("データが見つかりません。")


def chemical_space_smifp_comparison_display():
    st.title("📊 異なる化学空間の比較 (smifp)")

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

    # SELFIES対応表データ
    df_selfies = pd.DataFrame(data)

    st.write("ここでは、SMILESをカウントする代わりにSLFiesで変換してで比較し、図示します。")

    # Streamlitアプリの表示
    with st.expander("SMIfp項目とSELFIESトークンの対応表を表示"):
        st.write("以下は、SMIfp項目とSELFIESトークンの対応表です。")
        st.dataframe(df_selfies, use_container_width=True)



    # データの読み込み
    df = load_data()

    # selfiesのカウント
    try:
        from logic.selfies_features import count_selfies_features
        selfies_available = True
        st.success("SELFIES特徴量抽出ライブラリが利用可能です")
    except ImportError:
        selfies_available = False
        st.error("SELFIES特徴量抽出ライブラリが見つかりません")
        return

    if len(df) > 0 and selfies_available:
        # 処理を軽くするためサンプルデータに限定
        sample_df = df.sample(n=min(200, len(df)), random_state=42).copy()
        
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
            
            # 分子画像を生成（ホバー表示用）
            valid_sample_df['mol_image'] = valid_sample_df['Mol'].apply(mol_to_base64_image)
            
            # SELFIES特徴量の統計情報を表示
            st.subheader("📊 SELFIES特徴量の統計情報")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("処理した分子数", len(valid_sample_df))
                st.metric("特徴量次元数", len(selfies_df.columns))
            
            with col2:
                st.metric("TasteDB分子数", len(valid_sample_df[valid_sample_df['Source'] == 'TasteDB']))
                st.metric("FragranceDB分子数", len(valid_sample_df[valid_sample_df['Source'] == 'FragranceDB']))
            
            with col3:
                st.metric("平均分子量", f"{valid_sample_df['MW'].mean():.1f} Da")
                st.metric("平均LogP", f"{valid_sample_df['LogP'].mean():.2f}")
            
            # SELFIES特徴量の詳細表示
            with st.expander("🔍 SELFIES特徴量の詳細"):
                st.write("最初の10分子のSELFIES特徴量:")
                display_selfies = selfies_df.head(10).round(2)
                st.dataframe(display_selfies, use_container_width=True)
                
                # 特徴量の説明
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
                """)
            
            # PCAで2次元に次元削減
            features_for_pca = selfies_df.select_dtypes(include=[np.number]).fillna(0)
            
            if features_for_pca.shape[1] >= 2:
                pca = PCA(n_components=2)
                selfies_pca = pca.fit_transform(features_for_pca)
                
                valid_sample_df['SELFIES_PC1'] = selfies_pca[:, 0]
                valid_sample_df['SELFIES_PC2'] = selfies_pca[:, 1]
                
                # PCA結果の可視化
                st.subheader("📈 SELFIES特徴量空間での化学空間比較（PCA）")
                
                fig = go.Figure()
                
                sources = valid_sample_df['Source'].unique()
                colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
                
                for i, source in enumerate(sources):
                    source_data = valid_sample_df[valid_sample_df['Source'] == source]
                    
                    hover_template = (
                        "<b>%{customdata[0]}</b><br>"
                        "SMILES: %{customdata[1]}<br>"
                        "SELFIES PC1: %{x:.2f}<br>"
                        "SELFIES PC2: %{y:.2f}<br>"
                        "<extra></extra>"
                    )
                    
                    fig.add_trace(go.Scatter(
                        x=source_data['SELFIES_PC1'],
                        y=source_data['SELFIES_PC2'],
                        mode='markers',
                        name=source,
                        marker=dict(
                            color=colors[i % len(colors)],
                            size=8,
                            opacity=0.7
                        ),
                        customdata=source_data[['Source', 'SMILES', ]].values,
                        hovertemplate=hover_template
                    ))
                
                fig.update_layout(
                    title=f'SELFIES特徴量化学空間比較（PCA）- 寄与率: PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}',
                    xaxis_title=f'SELFIES PC1 ({pca.explained_variance_ratio_[0]:.1%})',
                    yaxis_title=f'SELFIES PC2 ({pca.explained_variance_ratio_[1]:.1%})',
                    width=800,
                    height=600,
                    hovermode='closest'
                )
                
                # Streamlitでプロット表示
                selected_points = st.plotly_chart(fig, use_container_width=True, on_select="rerun")
                
                # PCA寄与率の可視化
                st.subheader("📊 PCA寄与率の分析")
                
                # 各主成分の寄与率を計算（最大10成分まで）
                n_components_to_show = min(10, features_for_pca.shape[1])
                pca_full = PCA(n_components=n_components_to_show)
                pca_full.fit(features_for_pca)
                
                # 寄与率のバーチャート
                col1, col2 = st.columns(2)
                
                with col1:
                    # 個別寄与率
                    fig_contrib = go.Figure()
                    fig_contrib.add_trace(go.Bar(
                        x=[f'PC{i+1}' for i in range(n_components_to_show)],
                        y=pca_full.explained_variance_ratio_ * 100,
                        marker_color='lightblue',
                        text=[f'{ratio:.1f}%' for ratio in pca_full.explained_variance_ratio_ * 100],
                        textposition='auto'
                    ))
                    fig_contrib.update_layout(
                        title='各主成分の寄与率',
                        xaxis_title='主成分',
                        yaxis_title='寄与率 (%)',
                        height=400
                    )
                    st.plotly_chart(fig_contrib, use_container_width=True)
                
                with col2:
                    # 累積寄与率
                    cumulative_ratio = np.cumsum(pca_full.explained_variance_ratio_) * 100
                    fig_cumulative = go.Figure()
                    fig_cumulative.add_trace(go.Scatter(
                        x=[f'PC{i+1}' for i in range(n_components_to_show)],
                        y=cumulative_ratio,
                        mode='lines+markers',
                        marker_color='orange',
                        line=dict(width=3)
                    ))
                    fig_cumulative.add_hline(y=80, line_dash="dash", line_color="red", 
                                           annotation_text="80%ライン", annotation_position="bottom right")
                    fig_cumulative.update_layout(
                        title='累積寄与率',
                        xaxis_title='主成分',
                        yaxis_title='累積寄与率 (%)',
                        height=400
                    )
                    st.plotly_chart(fig_cumulative, use_container_width=True)
                
                # 寄与率の統計情報
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("PC1寄与率", f"{pca.explained_variance_ratio_[0]:.1%}")
                with col2:
                    st.metric("PC2寄与率", f"{pca.explained_variance_ratio_[1]:.1%}")
                with col3:
                    st.metric("PC1+PC2累積寄与率", f"{pca.explained_variance_ratio_[:2].sum():.1%}")

                # 主成分負荷量（特徴量の重み）を取得
                feature_names = features_for_pca.columns.tolist()
                pc1_loadings = pca.components_[0]
                pc2_loadings = pca.components_[1]
                
                # 寄与度の絶対値でソート
                pc1_importance = pd.DataFrame({
                    'feature': feature_names,
                    'loading': pc1_loadings,
                    'abs_loading': np.abs(pc1_loadings)
                }).sort_values('abs_loading', ascending=False)
                
                pc2_importance = pd.DataFrame({
                    'feature': feature_names,
                    'loading': pc2_loadings,
                    'abs_loading': np.abs(pc2_loadings)
                }).sort_values('abs_loading', ascending=False)

                # 寄与度上位を表示


                # 選択された分子の詳細情報
                st.subheader("🔍 選択された分子の詳細情報")
                
                # データソースフィルタリング
                source_filter = st.selectbox(
                    "データソースでフィルタ:",
                    ["全て"] + list(valid_sample_df['Source'].unique()),
                    key="selfies_source_filter"
                )
                
                # フィルタリングを適用
                filtered_df = valid_sample_df.copy()
                
                if source_filter != "全て":
                    filtered_df = filtered_df[filtered_df['Source'] == source_filter]
                
                st.write(f"フィルタ後の分子数: {len(filtered_df)} / {len(valid_sample_df)}")
                
                if len(filtered_df) > 0:
                    # 分子を選択するためのセレクトボックス
                    available_ids = filtered_df.index.tolist()
                    
                    selected_id = st.selectbox(
                        f"分子IDを選択して構造を確認 ({len(available_ids)} 件表示):", 
                        available_ids,
                        format_func=lambda x: f"ID: {x} - {filtered_df.loc[x, 'Source']} - MW: {filtered_df.loc[x, 'MW']:.1f} - SMILES: {filtered_df.loc[x, 'SMILES'][:30]}...",
                        key="selfies_molecule_selector"
                    )
                    
                    # 選択された分子の情報を表示
                    selected_mol = filtered_df.loc[selected_id]
                    selected_selfies_features = selfies_df.loc[selected_id]
                    
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
                        st.write(f"- PC1: {selected_mol['SELFIES_PC1']:.3f}")
                        st.write(f"- PC2: {selected_mol['SELFIES_PC2']:.3f}")
                    
                    with col2:
                        st.write("**分子構造:**")
                        if selected_mol['Mol'] is not None:
                            # RDKitで直接画像を生成してStreamlitで表示
                            mol_img = Draw.MolToImage(selected_mol['Mol'], size=(300, 300))
                            st.image(mol_img, caption="分子構造", use_container_width=True)
                        else:
                            st.error("分子構造を生成できませんでした")
                    
                    # SELFIES特徴量の詳細表示
                    st.write("**SELFIES特徴量:**")
                    selfies_col1, selfies_col2 = st.columns(2)
                    
                    with selfies_col1:
                        st.write("原子数:")
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
                        st.write("構造特徴:")
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
                
                else:
                    st.warning("フィルタ条件に一致する分子が見つかりません。")
                
                st.info("💡 SELFIES特徴量を使用することで、分子構造の特徴をより詳細に分析できます。PCAプロットとマウスオーバーで詳細情報を確認してください。")
                
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