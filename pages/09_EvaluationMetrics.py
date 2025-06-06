import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義

import streamlit as st
import pandas as pd
from logic.metrics import calculate_regression_metrics, create_yyplot
from logic.metrics import calculate_classification_metrics, calculate_confusion_matrix
from sklearn.datasets import load_diabetes

def split_data_display():
    """
    scikit-learnのtrain_test_splitまたはcross_val_scoreを使い、
    Streamlit上でホールドアウト分割か交差検証を選択し、分割割合やfold数を調整しながら
    統計情報やスコアを表示します。
    """
    import matplotlib.pyplot as plt
    from sklearn.datasets import load_iris
    from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
    from sklearn.ensemble import RandomForestClassifier

    st.title("📊 データ分割方法の解説")
    st.markdown("scikit-learnの`train_test_split`または`cross_val_score`を使ってデータセットを分割します。")

    # 分割方法の選択
    split_method = st.radio("分割方法を選択してください", ("Holdout", "Cross Validation"))

    # データセット選択
    st.subheader("データセット選択")
    dataset_source = st.radio("データソースを選択してください", ("Iris", "CSVファイルをアップロード"))
    if dataset_source == "Iris":
        from sklearn.datasets import load_iris
        data = load_iris(as_frame=True)
        X = data.data
        y = data.target

        # Irisデータセットの説明変数とtargetの値を表示
        st.subheader("Irisデータセット（特徴量とターゲット）")
        iris_df = X.copy()
        iris_df["target"] = y
        st.write(iris_df)
    else:
        uploaded_file = st.file_uploader("CSVファイルをアップロードしてください", type=["csv"])
        X, y = None, None
        if uploaded_file is not None:
            df = pd.read_csv(uploaded_file)
            st.write("アップロードされたデータ:")
            st.write(df)
            # target列を指定
            columns = list(df.columns)
            target_col = st.selectbox("目的変数（target）となる列を選択してください", columns)
            feature_cols = [col for col in columns if col != target_col]
            X = df[feature_cols]
            y = df[target_col]
            st.subheader("特徴量とターゲット")
            show_df = X.copy()
            show_df["target"] = y
            st.write(show_df)
        else:
            st.info("CSVファイルをアップロードしてください。")
            st.stop()

    if X is None or y is None:
        st.warning("データが正しく読み込まれていません。")
        st.stop()

    if split_method == "Holdout":
        # 分割割合と乱数シードの設定
        st.subheader("分割割合と乱数シードの設定")
        test_size = st.slider("Test size", min_value=0.05, max_value=0.5, value=0.2, step=0.05)
        val_size = st.slider("Validation size", min_value=0.05, max_value=0.5, value=0.2, step=0.05)
        random_seed = st.number_input("Random seed", min_value=0, max_value=10000, value=42, step=1)

        # まずテストデータを分割
        X_trainval, X_test, y_trainval, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_seed, stratify=y
        )
        # 残りから検証データを分割
        val_relative_size = val_size / (1 - test_size)
        X_train, X_val, y_train, y_val = train_test_split(
            X_trainval, y_trainval, test_size=val_relative_size, random_state=random_seed, stratify=y_trainval
        )

        # 統計情報の表示
        st.subheader("Dataset Statistics")
        data = {
            "Dataset": ["Train", "Validation", "Test"],
            "Rows": [len(X_train), len(X_val), len(X_test)],
            "Class 0": [
                sum(y_train == 0), sum(y_val == 0), sum(y_test == 0)
            ],
            "Class 1": [
                sum(y_train == 1), sum(y_val == 1), sum(y_test == 1)
            ],
            "Class 2": [
                sum(y_train == 2), sum(y_val == 2), sum(y_test == 2)
            ],
        }
        df = pd.DataFrame(data)
        st.table(df)

        # データインデックスの割り当て表示
        st.subheader("Sample Index Assignment")
        # 特徴量・目的変数・割り当てをまとめて表示
        index_df = X.copy()
        index_df["target"] = y
        index_df["Assignment"] = [
            "Train" if i in X_train.index else
            "Validation" if i in X_val.index else
            "Test" if i in X_test.index else
            "Unknown"
            for i in X.index
        ]
        st.write(index_df)

        # 円グラフ
        fig2, ax2 = plt.subplots()
        ax2.pie(df["Rows"], labels=df["Dataset"], autopct='%1.1f%%', colors=['blue', 'orange', 'green'])
        ax2.set_title("Dataset Proportion")
        st.pyplot(fig2)

    else:  # Cross Validation
        st.subheader("Cross Validation Settings")
        test_size = st.slider("Test size", min_value=0.05, max_value=0.5, value=0.2, step=0.05)
        random_seed = st.number_input("Random seed (CV)", min_value=0, max_value=10000, value=42, step=1)
        n_splits = st.slider("Number of folds (K)", min_value=2, max_value=10, value=5, step=1)

        # まずテストデータを分割
        from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score
        from sklearn.ensemble import RandomForestClassifier

        X_trainval, X_test, y_trainval, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_seed, stratify=y
        )

        # Cross Validation
        clf = RandomForestClassifier(random_state=42)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_seed)
        scores = cross_val_score(clf, X_trainval, y_trainval, cv=cv)

        st.write(f"Cross-validation scores for each fold: {scores}")
        st.write(f"Mean score: {scores.mean():.4f}")
        st.write(f"Std: {scores.std():.4f}")

        # 各サンプルが各foldで訓練/検証どちらに使われたかを記録
        assignment = pd.DataFrame(index=X_trainval.index)
        for fold, (train_idx, val_idx) in enumerate(cv.split(X_trainval, y_trainval), 1):
            assignment[f"Fold_{fold}"] = "Train"
            assignment.iloc[val_idx, assignment.columns.get_loc(f"Fold_{fold}")] = "Validation"
        # テストデータは常に"Test"
        assignment_full = pd.DataFrame(index=X.index)
        assignment_full["Assignment"] = [
            "Test" if i in X_test.index else "Train/Validation" for i in X.index
        ]
        assignment_full = assignment_full.join(assignment, how="left")

        st.subheader("各サンプルのFoldごとの割り当て")
        st.write(assignment_full)

        # 円グラフの作成
        st.subheader("データセット割合（円グラフ）")
        trainval_count = len(X_trainval)
        test_count = len(X_test)
        fig_cv, ax_cv = plt.subplots()
        ax_cv.pie(
            [trainval_count, test_count],
            labels=["Train/Validation", "Test"],
            autopct='%1.1f%%',
            colors=['blue', 'green']
        )
        ax_cv.set_title("Dataset Proportion (CV)")
        st.pyplot(fig_cv)





def display_regression_metrix():
    # ヘッダー
    st.title("📊 回帰評価指標計算ツール")
    st.markdown("アップロードされたCSVファイルやscikit-learnのデータセットを使用して回帰評価指標を計算します。")

    # データ入力セクション
    st.header("データ入力")
    data_option = st.radio(
        "データソースを選択してください：",
        ("scikit-learnのサンプルデータ", "CSVファイルをアップロード"),
        key="regression_data_source_radio"
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
            # CSVが無効な場合はreturn
            if not ({"actual", "predicted"}.issubset(df.columns)):
                return
        else:
            return
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

def display_classification_metrix():
    st.title("📊 分類評価指標計算ツール")
    st.markdown("アップロードされたCSVファイルやscikit-learnのデータセットを使用して分類評価指標を計算します。")

    # データ入力セクション
    st.header("データ入力")
    data_option = st.radio(
        "データソースを選択してください：",
        ("scikit-learnのサンプルデータ", "CSVファイルをアップロード"),
        key="classification_data_source_radio"
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
            if not ({"actual", "predicted"}.issubset(df.columns)):
                return
        else:
            return
    else:
        dataset_name = st.selectbox("Select a scikit-learn dataset", ["Iris (2 classes: 0/1)", ])
        if dataset_name == "Iris (2 classes: 0/1)":
            from sklearn.datasets import load_iris
            import numpy as np
            # Streamlitでrandom seedを指定
            random_seed = st.number_input("Random seed for prediction", min_value=0, max_value=10000, value=42, step=1)
            data = load_iris()
            # Use only class 0 and 1
            mask = data.target < 2
            y_true = data.target[mask]
            # Generate random predictions (0 or 1) with user seed
            rng = np.random.default_rng(seed=int(random_seed))
            y_pred = rng.integers(low=0, high=2, size=len(y_true))
            df = pd.DataFrame({"actual": y_true, "predicted": y_pred})
            st.write("Dataset used (Iris 2 classes):")
            st.write(df)

    # Calculate metrics
    if "y_true" in locals() and "y_pred" in locals():
        results = calculate_classification_metrics(y_true, y_pred)
        cm = calculate_confusion_matrix(y_true, y_pred)

        st.header("📈 Results")
        for metric, value in results.items():
            st.write(f"**{metric}**: {value:.4f}")

        st.header("📊 Confusion Matrix")
        labels = ["0", "1"]
        cm_df = pd.DataFrame(cm, index=[f"Actual:{l}" for l in labels], columns=[f"Predicted:{l}" for l in labels])
        st.dataframe(cm_df)

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "EvaluationMetrics"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()