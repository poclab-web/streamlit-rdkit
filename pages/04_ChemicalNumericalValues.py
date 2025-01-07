import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義
from logic.load_data import load_dataset_by_name

import streamlit as st
import seaborn as sns
import pandas as pd
from sklearn.datasets import load_iris, load_wine


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



if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ChemicalNumericalValues"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()