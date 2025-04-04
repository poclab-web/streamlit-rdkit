import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic.chemical_search import search_exact_match

import random
import pandas as pd

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

    data_file = 'data/Reagents/TCI_output_part1.csv'  # 修正済みのパス
    tci_data = pd.read_csv(data_file)
    smiles_list = tci_data['SMILES'].tolist()

    # ユーザーの入力
    query_smiles = st.text_input("検索したいSMILESを入力してください", "")

    # オプション設定
    ignore_stereo = st.checkbox("立体異性体を無視する", value=False)
    include_salts = st.checkbox("塩を含める（標準化しない）", value=True)

    if st.button("検索を実行"):
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
    st.title('RDKit + Py3DMOL 😀')

    # smartsを入力
    st.text("TCIで売られている化合物(数万個)の化合物をSMARTSで検索")
    search_smarts = st.text_input('SMARTSを入力', 'c1cc([Oh])ccc1')
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

    molsPerRow = st.text_input('構造例として１行に表示させる個数', '3')
    subImgSize = (300, 200)
    number = st.text_input('構造例として表示させたい分子数', '6')

    random_matches = random.sample(matches, int(number))
    img2 = Draw.MolsToGridImage(random_matches, molsPerRow=int(molsPerRow))
    st.image(img2)

    smi2 = [Chem.MolToSmiles(mol) for mol in random_matches]
    df2 = pd.DataFrame(smi2, columns=["smiles"])
    st.dataframe(df2)

    if st.button('合致した構造全てのcsvファイルの作成'):
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



if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "StructureSearch"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()