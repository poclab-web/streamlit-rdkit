import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from sklearn.datasets import load_wine
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, accuracy_score #回帰用の評価指標
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, ConfusionMatrixDisplay # 分類用の評価指標
from sklearn.preprocessing import StandardScaler

from sklearn.linear_model import LinearRegression, LogisticRegression # 線形回帰とロジスティック回帰
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier # K近傍法（回帰と分類）
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier # 決定木（回帰と分類）

import matplotlib.pyplot as plt
import pandas as pd

# データ取得用関数
def get_wine_data(task="regression", standardize=False):
    """
    task: "regression" または "classification"
    standardize: TrueならStandardScalerで標準化
    regression: color_intensityを目的変数とし、color_intensityとod280/od315を除く特徴量を使う
    classification: クラス0と1のみを使用、全特徴量使用
    """
    data = load_wine()
    X = pd.DataFrame(data.data, columns=data.feature_names)

    if task == "regression":
        y = X["color_intensity"]
        X = X.drop(["color_intensity", "od280/od315_of_diluted_wines", "hue"], axis=1)

    elif task == "classification":
        y = data.target
        mask = (y == 0) | (y == 1)
        X = X[mask]
        y = y[mask]
    else:
        raise ValueError("taskは'regression'または'classification'で指定してください。")

    if standardize:
        scaler = StandardScaler()
        X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)

    return X, y

def linear_regression_display():
    st.title("線形回帰")
    st.write("ワインデータセットの、color intensityを予測します。")

    # ⭐ 標準化オプションのチェックボックス
    use_standardize = st.checkbox("特徴量を標準化する", value=True)

    # データ取得（標準化オプション反映）
    X, y = get_wine_data(task="regression", standardize=use_standardize)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    model = LinearRegression()
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標の計算
    mse_train = mean_squared_error(y_train, y_pred_train)
    rmse_train = mse_train ** 0.5
    mae_train = (abs(y_train - y_pred_train)).mean()
    r2_train = model.score(X_train, y_train)

    mse_test = mean_squared_error(y_test, y_pred_test)
    rmse_test = mse_test ** 0.5
    mae_test = (abs(y_test - y_pred_test)).mean()
    r2_test = model.score(X_test, y_test)

    metrics_df = pd.DataFrame({
        "RMSE": [rmse_train, rmse_test],
        "MAE": [mae_train, mae_test],
        "R2": [r2_train, r2_test]
    }, index=["Train", "Test"])
    st.subheader("評価指標")
    st.table(metrics_df)

    # yyプロット
    
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_pred_train, alpha=0.7, label="Train", color="blue")
    ax.scatter(y_test, y_pred_test, alpha=0.7, label="Test", color="orange")
    min_y = min(y.min(), y_pred_train.min(), y_pred_test.min())
    max_y = max(y.max(), y_pred_train.max(), y_pred_test.max())
    ax.plot([min_y, max_y], [min_y, max_y], 'r--')
    ax.set_xlabel("True color intensity")
    ax.set_ylabel("Predicted color intensity")
    ax.set_title("yyplot")
    ax.legend()
    st.pyplot(fig)

    st.subheader("回帰式")
    coef_df = pd.DataFrame({
        "特徴量": X.columns,
        "係数": model.coef_
    })
    st.table(coef_df)
    st.write(f"切片（intercept）: {model.intercept_:.3f}")


def logistic_regression_display():
    st.title("ロジスティック回帰")
    st.write("ワインデータセットのうちクラス0と1を使って、ロジスティック回帰で分類を行います。")

    # 標準化のチェックボックス
    use_standardize = st.checkbox("特徴量を標準化する", value=True, key="logreg_std")

    # データ取得
    X, y = get_wine_data(task="classification", standardize=use_standardize)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # モデル学習
    model = LogisticRegression(random_state=0)
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標の計算
    metrics_df = pd.DataFrame({
        "Accuracy": [accuracy_score(y_train, y_pred_train), accuracy_score(y_test, y_pred_test)],
        "Precision": [precision_score(y_train, y_pred_train), precision_score(y_test, y_pred_test)],
        "Recall": [recall_score(y_train, y_pred_train), recall_score(y_test, y_pred_test)],
        "F1 Score": [f1_score(y_train, y_pred_train), f1_score(y_test, y_pred_test)],
    }, index=["Train", "Test"])

    st.subheader("評価指標")
    st.table(metrics_df)

    # 混同行列（TrainとTestの両方）
    st.subheader("混同行列")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Trainデータ**")
        cm_train = confusion_matrix(y_train, y_pred_train)
        fig1, ax1 = plt.subplots()
        disp1 = ConfusionMatrixDisplay(confusion_matrix=cm_train, display_labels=model.classes_)
        disp1.plot(ax=ax1, cmap="Blues", colorbar=False)
        st.pyplot(fig1)

    with col2:
        st.markdown("**Testデータ**")
        cm_test = confusion_matrix(y_test, y_pred_test)
        fig2, ax2 = plt.subplots()
        disp2 = ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=model.classes_)
        disp2.plot(ax=ax2, cmap="Oranges", colorbar=False)
        st.pyplot(fig2)

    # 回帰係数の表示
    st.subheader("回帰係数（特徴量ごとの重み）")
    coef_df = pd.DataFrame({
        "特徴量": X.columns,
        "係数": model.coef_[0]
    })
    st.table(coef_df)
    st.write(f"切片（intercept）: {model.intercept_[0]:.3f}")


def knn_regression_display():
    st.title("KNN回帰")
    st.write("ワインデータセットのデータを使い、color intensityをK近傍法（KNN）で予測します。")

    # 標準化チェックボックス
    use_standardize = st.checkbox("特徴量を標準化する", value=True, key="knn_std")

    # 最近傍の数を選択（ハイパーパラメータ）
    k = st.slider("近傍数 k", min_value=1, max_value=20, value=5, step=1)

    # データ取得
    X, y = get_wine_data(task="regression", standardize=use_standardize)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # モデル構築・学習
    model = KNeighborsRegressor(n_neighbors=k)
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標の計算
    mse_train = mean_squared_error(y_train, y_pred_train)
    rmse_train = mse_train ** 0.5
    mae_train = (abs(y_train - y_pred_train)).mean()
    r2_train = model.score(X_train, y_train)

    mse_test = mean_squared_error(y_test, y_pred_test)
    rmse_test = mse_test ** 0.5
    mae_test = (abs(y_test - y_pred_test)).mean()
    r2_test = model.score(X_test, y_test)

    metrics_df = pd.DataFrame({
        "RMSE": [rmse_train, rmse_test],
        "MAE": [mae_train, mae_test],
        "R2": [r2_train, r2_test]
    }, index=["Train", "Test"])
    st.subheader("評価指標")
    st.table(metrics_df)

    # yyプロット
    
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_pred_train, alpha=0.7, label="Train", color="blue")
    ax.scatter(y_test, y_pred_test, alpha=0.7, label="Test", color="orange")
    min_y = min(y.min(), y_pred_train.min(), y_pred_test.min())
    max_y = max(y.max(), y_pred_train.max(), y_pred_test.max())
    ax.plot([min_y, max_y], [min_y, max_y], 'r--')
    ax.set_xlabel("True color intensity")
    ax.set_ylabel("Predicted color intensity")
    ax.set_title(f"yyplot (k={k})")
    ax.legend()
    st.pyplot(fig)

def knn_classification_display():
    st.title("KNN分類")
    st.write("ワインデータセットのうちクラス0と1を使って、K近傍法（KNN）で分類を行います。")

    # 標準化のチェックボックス（key付きで重複防止）
    use_standardize = st.checkbox("特徴量を標準化する", value=True, key="knn_class_std")

    # k（近傍数）のスライダー（key付きで重複防止）
    k = st.slider("近傍数 k", min_value=1, max_value=20, value=5, step=1, key="knn_class_k")

    # データ取得
    X, y = get_wine_data(task="classification", standardize=use_standardize)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # モデル学習
    model = KNeighborsClassifier(n_neighbors=k)
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標の計算
    metrics_df = pd.DataFrame({
        "Accuracy": [accuracy_score(y_train, y_pred_train), accuracy_score(y_test, y_pred_test)],
        "Precision": [precision_score(y_train, y_pred_train), precision_score(y_test, y_pred_test)],
        "Recall": [recall_score(y_train, y_pred_train), recall_score(y_test, y_pred_test)],
        "F1 Score": [f1_score(y_train, y_pred_train), f1_score(y_test, y_pred_test)],
    }, index=["Train", "Test"])

    st.subheader("評価指標")
    st.table(metrics_df)

    # 混同行列（TrainとTestの両方）
    st.subheader("混同行列")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Trainデータ**")
        cm_train = confusion_matrix(y_train, y_pred_train)
        fig1, ax1 = plt.subplots()
        disp1 = ConfusionMatrixDisplay(confusion_matrix=cm_train, display_labels=model.classes_)
        disp1.plot(ax=ax1, cmap="Blues", colorbar=False)
        st.pyplot(fig1)

    with col2:
        st.markdown("**Testデータ**")
        cm_test = confusion_matrix(y_test, y_pred_test)
        fig2, ax2 = plt.subplots()
        disp2 = ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=model.classes_)
        disp2.plot(ax=ax2, cmap="Oranges", colorbar=False)
        st.pyplot(fig2)


def decision_tree_regression_display():
    st.title("決定木回帰")
    st.write("ワインデータセットの全データを使い、color intensityを決定木で予測します。")

    # 特徴量の標準化オプション（決定木では通常不要だが比較のため）
    use_standardize = st.checkbox("特徴量を標準化する", value=False, key="tree_std")

    # max_depth を選べるようにする
    max_depth = st.slider("最大深さ (max_depth)", min_value=1, max_value=20, value=5, step=1)

    # データの取得
    X, y = get_wine_data(task="regression", standardize=use_standardize)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # 決定木モデルの構築と学習
    model = DecisionTreeRegressor(max_depth=max_depth, random_state=0)
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標の計算
    mse_train = mean_squared_error(y_train, y_pred_train)
    rmse_train = mse_train ** 0.5
    mae_train = (abs(y_train - y_pred_train)).mean()
    r2_train = model.score(X_train, y_train)

    mse_test = mean_squared_error(y_test, y_pred_test)
    rmse_test = mse_test ** 0.5
    mae_test = (abs(y_test - y_pred_test)).mean()
    r2_test = model.score(X_test, y_test)

    # 表形式で表示
    metrics_df = pd.DataFrame({
        "RMSE": [rmse_train, rmse_test],
        "MAE": [mae_train, mae_test],
        "R2": [r2_train, r2_test]
    }, index=["Train", "Test"])
    st.subheader("評価指標")
    st.table(metrics_df)

    # yyプロット
    
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_pred_train, alpha=0.7, label="Train", color="blue")
    ax.scatter(y_test, y_pred_test, alpha=0.7, label="Test", color="orange")
    min_y = min(y.min(), y_pred_train.min(), y_pred_test.min())
    max_y = max(y.max(), y_pred_train.max(), y_pred_test.max())
    ax.plot([min_y, max_y], [min_y, max_y], 'r--')
    ax.set_xlabel("True color intensity")
    ax.set_ylabel("Predicted color intensity")
    ax.set_title(f"yyplot (max_depth={max_depth})")
    ax.legend()
    st.pyplot(fig)


def decision_tree_classification_display():
    st.title("決定木分類")
    st.write("ワインデータセットのうちクラス0と1を使って、決定木で分類を行います。")

    # 標準化チェックボックス（通常不要だが比較のため追加）
    use_standardize = st.checkbox("特徴量を標準化する", value=False, key="tree_class_std")

    # max_depthスライダー
    max_depth = st.slider("木の最大深さ (max_depth)", min_value=1, max_value=20, value=5, step=1, key="tree_class_depth")

    # データ取得
    X, y = get_wine_data(task="classification", standardize=use_standardize)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # モデル学習
    model = DecisionTreeClassifier(max_depth=max_depth, random_state=0)
    model.fit(X_train, y_train)
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)

    # 評価指標
    metrics_df = pd.DataFrame({
        "Accuracy": [accuracy_score(y_train, y_pred_train), accuracy_score(y_test, y_pred_test)],
        "Precision": [precision_score(y_train, y_pred_train), precision_score(y_test, y_pred_test)],
        "Recall": [recall_score(y_train, y_pred_train), recall_score(y_test, y_pred_test)],
        "F1 Score": [f1_score(y_train, y_pred_train), f1_score(y_test, y_pred_test)],
    }, index=["Train", "Test"])

    st.subheader("評価指標")
    st.table(metrics_df)

    # 混同行列
    st.subheader("混同行列")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Trainデータ**")
        cm_train = confusion_matrix(y_train, y_pred_train)
        fig1, ax1 = plt.subplots()
        disp1 = ConfusionMatrixDisplay(confusion_matrix=cm_train, display_labels=model.classes_)
        disp1.plot(ax=ax1, cmap="Blues", colorbar=False)
        st.pyplot(fig1)

    with col2:
        st.markdown("**Testデータ**")
        cm_test = confusion_matrix(y_test, y_pred_test)
        fig2, ax2 = plt.subplots()
        disp2 = ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=model.classes_)
        disp2.plot(ax=ax2, cmap="Oranges", colorbar=False)
        st.pyplot(fig2)


if __name__ == "__main__":
    current_category = "MachineLearning"
    st.write(f"現在のカテゴリー: {current_category}")

    handle_tabs_for_category(current_category)
    display_sidebar()