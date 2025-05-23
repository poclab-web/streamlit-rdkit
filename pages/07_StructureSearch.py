import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic.chemical_search import search_exact_match
from logic.similarity import find_similar_compounds, calculate_similarity

import random
import pandas as pd
import concurrent.futures
import time

from rdkit import Chem
from rdkit.Chem import Draw

# アプリの定義

## 部分構造検索
@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')

def search_exact_match_display():
    # Streamlitアプリ
    st.title("化合物データ検索アプリ")
    
    st.text("TCIで売られている化合物(数万個)の化合物を構造一致検索")

    data_file = 'data/Reagents/TCI_output_part1.csv'  # 修正済みのパス
    tci_data = pd.read_csv(data_file)
    smiles_list = tci_data['SMILES'].tolist()

    # ユーザーの入力
    query_smiles = st.text_input("検索したいSMILESを入力してください", "C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O", key="exact_match_query")

    # オプション設定
    ignore_stereo = st.checkbox("立体異性体を無視する", value=True)
    include_salts = st.checkbox("塩を含めず検索", value=False)

    if st.button("検索を実行", key="exact_match_search"):
        if not query_smiles:
            st.warning("SMILESを入力してください！")
        else:
            matches = search_exact_match(query_smiles, smiles_list, ignore_stereo=ignore_stereo, include_salts=include_salts)

            if matches:
                st.success(f"以下の化合物が見つかりました (合計: {len(matches)}):")
                for match in matches:
                    st.write(f"ヒットしたSMILES: {match}")
                    
                    # 構造式を表示
                    mol = Chem.MolFromSmiles(match)
                    st.image(Draw.MolToImage(mol))
                    
                    # 詳細情報の表示
                    matched_data = tci_data[tci_data['SMILES'] == match]

                    # 必要な情報のみ選択
                    selected_columns = [
                        'PUBCHEM_SUBSTANCE_SYNONYM', 
                        'PUBCHEM_EXT_SUBSTANCE_URL', 
                        'PUBCHEM_EXT_DATASOURCE_REGID',
                        'PUBCHEM_CID_ASSOCIATIONS',
                        'SMILES'
                    ]
                    if all(col in matched_data.columns for col in selected_columns):
                        for _, row in matched_data.iterrows():
                            # 詳細情報を見やすく表示
                            st.write(f"**PubChemID**: {row['PUBCHEM_CID_ASSOCIATIONS']}")
                            st.write(f"**物質名**: {row['PUBCHEM_SUBSTANCE_SYNONYM']}")
                            st.write(f"**試薬番号**: {row['PUBCHEM_EXT_DATASOURCE_REGID']}")
                            st.markdown(
                                f"**URL**: [{row['PUBCHEM_EXT_SUBSTANCE_URL']}]({row['PUBCHEM_EXT_SUBSTANCE_URL']})", unsafe_allow_html=True
                            )
                            st.write(f"**SMILES**: {row['SMILES']}")
                    else:
                        st.warning("表示可能なカラムが見つかりませんでした。元のデータを確認してください。")
            else:
                st.error("該当する化合物は見つかりませんでした。")

def smarts_search_display():
    st.title('Smarts Search 😀')

    # smartsを入力
    st.text("TCIで売られている化合物(数万個)の化合物をSMARTSで検索")
    search_smarts = st.text_input('SMARTSを入力', 'c1cc([Oh])ccc1', key="smarts_query")
    patt = Chem.MolFromSmarts(search_smarts)

    # imgを表示
    st.text("smartsで読み込んだ画像表示")
    img = Draw.MolsToGridImage([patt])
    st.image(img)
    st.code(search_smarts)

    # TCIデータの読み込み
    df = pd.read_csv("data/TCI_smiles.csv", encoding='shift_jis', index_col=0)

    mols = []
    for smi in df["smiles"]:
        try:
            mol = Chem.MolFromSmiles(smi)
            mols.append(mol)
        except:
            pass

    matches = [mol for mol in mols if mol.HasSubstructMatch(patt)]
    st.markdown("### 合致した構造の数" + str(len(matches)))

    molsPerRow = st.text_input('構造例として１行に表示させる個数', '3', key="mols_per_row")
    subImgSize = (300, 200)
    number = st.text_input('構造例として表示させたい分子数', '6', key="number_of_molecules")

    random_matches = random.sample(matches, int(number))
    img2 = Draw.MolsToGridImage(random_matches, molsPerRow=int(molsPerRow))
    st.image(img2)

    smi2 = [Chem.MolToSmiles(mol) for mol in random_matches]
    df2 = pd.DataFrame(smi2, columns=["smiles"])
    st.dataframe(df2)

    if st.button('合致した構造全てのcsvファイルの作成', key="create_csv"):
        smiles_list = [Chem.MolToSmiles(mol) for mol in matches]
        df_download = pd.DataFrame(smiles_list, columns=["smiles"])
        df_download_csv = convert_df(df_download)

        st.download_button(
            "smiles.csvのDownload",
            df_download_csv,
            "smiles.csv",
            "text/csv",
            key='download-csv'
        )
    else:
        st.write('Please Click Start Download button!')


def similarity_search_display():
    st.title('Similarity Search 😀')

    # 注意書きを追加
    st.warning("※ 類似性検索では現在、立体異性体や塩の扱いについて特別な処理は行っていません。今後、処理を加えます。")

    data_file = 'data/Reagents/TCI_output_part1.csv'  
    tci_data = pd.read_csv(data_file)
    smiles_list = tci_data['SMILES'].tolist()

    # ユーザーの入力
    query_smiles = st.text_input("検索したいSMILESを入力してください", "C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O", key="similarity_query")

    # 表示件数
    top_n = st.number_input("類似性が近いものから何件表示しますか？", min_value=1, max_value=100, value=10, step=1)

    # 優先する手法を選択
    method_options = {
        "Tanimoto(ECFP radius2 2048bit)": "fingerprint",
        "Normalized Levenshtein": "levenshtein",
        "Descriptor": "descriptor",
        "MCS": "mcs"
    }
    score_key_map = {
        "fingerprint": "fingerprint_similarity",
        "levenshtein": "normalized_levenshtein_distance",
        "descriptor": "descriptor_distance",
        "mcs": "mcs_similarity"
    }
    method_label = st.selectbox("何を優先して並べますか？", list(method_options.keys()))
    selected_method = method_options[method_label]
    selected_score_key = score_key_map[selected_method]

    if st.button("検索を実行", key="similarity_search"):
        if not query_smiles:
            st.warning("SMILESを入力してください！")
        else:
            all_methods = list(method_options.values())
            results = []

            def calc(smiles):
                scores = calculate_similarity(query_smiles, smiles, methods=all_methods)
                return (smiles, scores)

            with concurrent.futures.ThreadPoolExecutor() as executor:
                calc_results = list(executor.map(calc, smiles_list))

            # 並べ替え方向（Tanimoto, MCSは降順、他は昇順）
            reverse_sort = True if selected_method in ["fingerprint", "mcs"] else False
            calc_results.sort(key=lambda x: x[1][selected_score_key], reverse=reverse_sort)

            # 上位N件のみ
            top_results = calc_results[:top_n]

            if top_results:
                st.success(f"上位 {top_n} 件を表示します（{method_label}優先）")
                # データフレームを作成して表示
                columns = ["SMILES"] + [score_key_map[m] for m in all_methods]
                df_similar = pd.DataFrame(
                    [
                        [smiles] + [scores[score_key_map[m]] for m in all_methods]
                        for smiles, scores in top_results
                    ],
                    columns=columns
                )
                st.dataframe(df_similar)
            else:
                st.error("該当する化合物は見つかりませんでした。")

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "StructureSearch"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()