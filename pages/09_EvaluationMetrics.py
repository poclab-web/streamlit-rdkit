import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義

import streamlit as st
import pandas as pd
from logic.metrics import calculate_regression_metrics, create_yyplot
from sklearn.datasets import load_diabetes

def display_regression_metrix():
    # ヘッダー
    st.title("📊 回帰評価指標計算ツール")
    st.markdown("アップロードされたCSVファイルやscikit-learnのデータセットを使用して回帰評価指標を計算します。")

    # データ入力セクション
    st.header("データ入力")
    data_option = st.radio(
        "データソースを選択してください：",
        ("CSVファイルをアップロード", "scikit-learnのサンプルデータ")
    )

    # データ読み込み
    if data_option == "CSVファイルをアップロード":
        uploaded_file = st.file_uploader("CSVファイルをアップロードしてください", type=["csv"])
        if uploaded_file is not None:
            df = pd.read_csv(uploaded_file)
            st.write("アップロードされたデータ:")
            st.write(df)
            if {"actual", "predicted"}.issubset(df.columns):
                y_true = df["actual"]
                y_pred = df["predicted"]
            else:
                st.error("CSVファイルには 'actual' 列と 'predicted' 列が必要です。")
    else:
        dataset_name = st.selectbox("scikit-learnのデータセットを選択してください", ["Diabetes", ])
        if dataset_name == "Diabetes":
            data = load_diabetes()

        X, y = data.data, data.target
        y_true = y[:50]  # 最初の50サンプルを実際の値とする
        y_pred = y_true + (0.1 * y_true.std()) * (2 * (pd.Series(range(50)) % 2) - 1)  # ノイズを加えた予測値
        df = pd.DataFrame({"actual": y_true, "predicted": y_pred})
        st.write("使用するデータセット:")
        st.write(df)

    # 評価指標の計算
    if "y_true" in locals() and "y_pred" in locals():
        results = calculate_regression_metrics(y_true, y_pred)

        # プロット作成と表示
        st.header("📊 y-y プロット")
        fig = create_yyplot(y_true, y_pred)  # プロット作成
        st.pyplot(fig)  # プロットをStreamlitで表示

        st.header("📈 計算結果")
        for metric, value in results.items():
            st.write(f"**{metric}**: {value:.4f}")


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "EvaluationMetrics"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()