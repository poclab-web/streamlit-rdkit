import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義
from logic.load_data import load_dataset_by_name
from logic.mol_loader import MoleculeDataLoader

import streamlit as st
import seaborn as sns
import pandas as pd
from sklearn.datasets import load_iris, load_wine

from rdkit.Chem import Draw

# 有名なデータセットの読み込み

def display_load_data():
    # タイトル
    st.title("有名なデータセットを切り替えて表示")

    # データセットの選択
    dataset_name = st.selectbox(
        "データセットを選択してください",
        ["Titanic", "Iris", "Penguins", "Wine"]
    )

    # 選択したデータセットをロード
    df = load_dataset_by_name(dataset_name)

    # データセットの中身表示
    st.write(f"### {dataset_name} データセットの中身")
    st.dataframe(df.head())

    # 各カラムのデータ型表示
    st.write("### 各カラムのデータ型")
    st.text(df.dtypes)

    # データセットの基本統計情報
    st.write("### データの基本統計情報")
    st.write(df.describe())

# 溶解度、融点、molLogPなど
def dataset_viewer():
    """
    Streamlitアプリ: SMILESファイルやデータセットを統一的にロードして表示。
    """
    st.title("Molecule Data Viewer")

    # アップロードセクション
    st.header("データアップロード")
    uploaded_file = st.file_uploader("ファイルをアップロード (SMILES/CSV/Excel)", type=["smi", "csv", "txt", "xlsx"])

    if uploaded_file:
        try:
            # ファイル拡張子に応じた処理
            if uploaded_file.name.endswith(".smi") or uploaded_file.name.endswith(".txt"):
                data = MoleculeDataLoader.load_smiles_file(uploaded_file)
            elif uploaded_file.name.endswith(".csv") or uploaded_file.name.endswith(".tab"):
                data = MoleculeDataLoader.load_csv(uploaded_file, separator='\t')
            elif uploaded_file.name.endswith(".xlsx"):
                data = MoleculeDataLoader.load_excel(uploaded_file)
            else:
                st.error("対応していないファイル形式です。")
                return

            # データ表示
            st.write("### アップロードデータ")
            st.dataframe(data)

            # 分子構造の表示オプション
            if st.checkbox("分子構造を表示"):
                mols = data["Mol"].tolist()
                legends = data["Name"].tolist() if "Name" in data.columns else None
                img = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=4)
                st.image(img)

        except Exception as e:
            st.error(f"データの処理中にエラーが発生しました: {e}")

    # サンプルデータのセクション
    st.header("サンプルデータ")
    if st.button("Load Sample Data"):
        try:
            # データセットの選択
            dataset_name = st.selectbox(
                "データセットを選択してください",
                ["solubility", "NMR", "food", "qssr"]
            )

            # サンプルファイルのパス
            sample_dataset_file_solubility = 'data/curated-solubility-dataset.tab'
            sample_smiles_files_food = ['data/TasteDB.smi', 'data/FragranceDB.smi']
            sample_dataset_file_NMR = 'data/NMRshiftDB2_CHOonly_no_missing.xlsx'


        except Exception as e:
            st.error(f"サンプルデータのロード中にエラーが発生しました: {e}")



if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ChemicalNumericalValues"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()