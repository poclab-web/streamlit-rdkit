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
    st.title("Molecule Data Viewer")

    # サンプルデータのセクション
    st.header("サンプルデータ")
    
    # 選択したデータセットをセッション状態で保持
    if "selected_dataset" not in st.session_state:
        st.session_state["selected_dataset"] = "solubility"

    # データセットの選択
    dataset_name = st.selectbox(
        "データセットを選択してください",
        ["solubility", "NMR", "food", "qssr"],
        key="selected_dataset"
    )

    if st.button("Load Sample Data"):
        try:
            # サンプルデータのロード処理
            if dataset_name == "solubility":
                file_path = 'data/curated-solubility-dataset.tab'
                data = pd.read_csv(file_path, sep='\\t')
            elif dataset_name == "NMR":
                file_path = 'data/NMRshiftDB2_CHOonly_no_missing.xlsx'
                data = pd.read_excel(file_path)
            elif dataset_name == "food":
                file_paths = ['data/TasteDB.smi', 'data/FragranceDB.smi']
                data_frames = [MoleculeDataLoader.load_smiles_file(file) for file in file_paths]
                data = pd.concat(data_frames, ignore_index=True)
            elif dataset_name == "qssr":
                st.warning("QSSRデータセットのロード機能は未実装です。")
                return
            else:
                st.error("不明なデータセットが選択されました。")
                return

            # サンプルデータの表示
            st.write(f"### サンプルデータ: {dataset_name}")
            st.dataframe(data)

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