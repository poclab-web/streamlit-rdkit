import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import NearestNeighbors

import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler

from scipy.interpolate import griddata

from sklearn.datasets import load_wine
import matplotlib.pyplot as plt

import seaborn as sns
import plotly.express as px
import numpy as np

import matplotlib.patches as patches
from matplotlib.patches import Circle


import shap

# アプリの定義
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


def visualization_display():
    # UI表示
    st.title("🍷 ワインデータの主成分分析 (PCA)")

    # データ読み込みと前処理
    X, y = get_wine_data(task="classification", standardize=True)
    class_labels = pd.Series(y).map({0: "0", 1: "1"})

    # カラーマッピング
    color_map = {"0": "blue", "1": "red"}
    colors = class_labels.map(color_map)

    # 2次元PCA
    pca_2d = PCA(n_components=2)
    X_pca_2d = pca_2d.fit_transform(X)

    fig1 = px.scatter(
        x=X_pca_2d[:, 0], y=X_pca_2d[:, 1],
        color=class_labels,
        color_discrete_map=color_map,
        labels={"x": "PC1", "y": "PC2", "color": "Class"},
        title="2次元PCAプロット"
    )
    st.plotly_chart(fig1)

    # 3次元PCA
    pca_3d = PCA(n_components=3)
    X_pca_3d = pca_3d.fit_transform(X)

    fig2 = px.scatter_3d(
        x=X_pca_3d[:, 0], y=X_pca_3d[:, 1], z=X_pca_3d[:, 2],
        color=class_labels,
        color_discrete_map=color_map,
        size_max=4,
        size=[3] * len(X),  # サイズを小さく
        labels={"x": "PC1", "y": "PC2", "z": "PC3", "color": "Class"},
        title="3次元PCAプロット"
    )
    st.plotly_chart(fig2)

    # 累積寄与率プロット
    pca_full = PCA(n_components=X.shape[1])
    pca_full.fit(X)
    exp_var = pca_full.explained_variance_ratio_
    cum_var = np.cumsum(exp_var)

    fig3, ax = plt.subplots()
    ax.bar(range(1, len(exp_var) + 1), exp_var, alpha=0.5, label="Individual")
    ax.step(range(1, len(cum_var) + 1), cum_var, where="mid", label="Cumulative")
    ax.set_xlabel("Principal Component Index")
    ax.set_ylabel("Explained Variance Ratio")
    ax.set_title("Explained Variance Ratio and Cumulative Contribution")
    ax.legend(loc="best")
    st.pyplot(fig3)

    # 主成分負荷量のヒートマップ
    loadings = pd.DataFrame(
        pca_full.components_.T,
        columns=[f"PC{i+1}" for i in range(len(pca_full.components_))],
        index=X.columns
    )

    st.subheader("🔍 Principal Component Loadings")
    fig4, ax2 = plt.subplots(figsize=(10, 6))
    sns.heatmap(loadings.iloc[:, :5], annot=True, cmap="coolwarm", center=0, ax=ax2)
    ax2.set_title("Loadings of Top 5 Principal Components")
    ax2.set_ylabel("Features")
    ax2.set_xlabel("Principal Components")
    st.pyplot(fig4)

def AD_visualization_display():
    st.title("🍷 PCA + Applicability Domain (Class 0 only)")

    # 1. データ取得とPCA
    X, y = get_wine_data(task="classification", standardize=True)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)

    # 2. データフレーム化
    df_pca = pd.DataFrame(X_pca, columns=["PC1", "PC2"])
    df_pca["Class"] = y

    # 3. クラス0の範囲でADを定義
    df_class0 = df_pca[df_pca["Class"] == 0]
    x_min, x_max = df_class0["PC1"].min(), df_class0["PC1"].max()
    y_min, y_max = df_class0["PC2"].min(), df_class0["PC2"].max()

    # 4. プロット
    fig5, ax = plt.subplots()
    for cls, color in zip([0, 1], ["blue", "red"]):
        ax.scatter(df_pca[df_pca["Class"] == cls]["PC1"],
                df_pca[df_pca["Class"] == cls]["PC2"],
                label=f"Class {cls}", color=color, alpha=0.6)

    # クラス0に基づく適用範囲（AD）を緑枠で描画
    rect = patches.Rectangle(
        (x_min, y_min), x_max - x_min, y_max - y_min,
        linewidth=2, edgecolor='green', facecolor='none',
        label="AD: Class 0"
    )
    ax.add_patch(rect)

    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_title("Ranging-based AD defined by Class 0 only")
    ax.legend()
    st.pyplot(fig5)



    # --- Streamlit UI：パーセンタイル選択 ---
    st.header("Distance-based AD 設定")
    percentile = st.slider("しきい値のパーセンタイル（%）", min_value=50, max_value=100, value=95, step=1)

    # --- 中心（重心）を計算 ---
    center_x = df_class0["PC1"].mean()
    center_y = df_class0["PC2"].mean()

    # --- 各点との距離を計算 ---
    df_pca["distance_to_center"] = np.sqrt((df_pca["PC1"] - center_x)**2 + (df_pca["PC2"] - center_y)**2)

    # --- 適用範囲の半径を選択したパーセンタイルで決定 ---
    threshold = np.percentile(df_class0.apply(
        lambda row: np.sqrt((row["PC1"] - center_x)**2 + (row["PC2"] - center_y)**2), axis=1
    ), percentile)

    # --- 可視化 ---
    fig6, ax = plt.subplots()

    # 各クラスの点をプロット
    for cls, color in zip([0, 1], ["blue", "red"]):
        ax.scatter(df_pca[df_pca["Class"] == cls]["PC1"],
                df_pca[df_pca["Class"] == cls]["PC2"],
                label=f"Class {cls}", color=color, alpha=0.5)

    # 重心とAD円を描画
    ax.plot(center_x, center_y, marker="*", color="red", markersize=15, label="Class 0 Center")
    circle = Circle((center_x, center_y), threshold, color="green", fill=False, linewidth=2, linestyle="--", label=f"AD radius ({percentile}%)")
    ax.add_patch(circle)

    # 軸・ラベル
    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_title(f"Distance-based AD (Class 0, {percentile} percentile)")
    ax.legend()
    st.pyplot(fig6)


    # --- データ準備（df_class0, df_class1, df_pca は定義済とする） ---
    X_train = df_class0[["PC1", "PC2"]].values
    df_class1 = df_pca[df_pca["Class"] == 1]

    # --- Streamlit UI（距離しきい値スライダー） ---
    st.header("KNN-AD 設定")
    k = st.slider("近傍数 k", min_value=1, max_value=10, value=3, step=1)
    threshold = st.slider("距離しきい値 (ADの境界)", min_value=0.1, max_value=0.5, value=0.4, step=0.025)

    # --- KNN距離の計算 ---
    nn_model = NearestNeighbors(n_neighbors=k)
    nn_model.fit(X_train)
    distances, _ = nn_model.kneighbors(X_train)
    knn_radii = distances[:, -1]

    # --- グリッド作成と補間 ---
    grid_x, grid_y = np.mgrid[
        X_train[:, 0].min() - 1:X_train[:, 0].max() + 1:200j,
        X_train[:, 1].min() - 1:X_train[:, 1].max() + 1:200j
    ]
    grid_z = griddata(X_train, knn_radii, (grid_x, grid_y), method='linear')
    grid_z[np.isnan(grid_z)] = np.nanmedian(grid_z)

    # --- 可視化 ---
    fig7, ax = plt.subplots()

    # Class 0, 1 のプロット
    ax.scatter(df_class0["PC1"], df_class0["PC2"], color="blue", label="Class 0", alpha=0.6)
    ax.scatter(df_class1["PC1"], df_class1["PC2"], color="red", label="Class 1", alpha=0.6)

    # # 等高線（連続距離線）
    # contour = ax.contour(grid_x, grid_y, grid_z, levels=6, colors="green", linewidths=1)
    # ax.clabel(contour, inline=True, fontsize=8)

    # 塗り分け：AD内（距離≦閾値）領域を緑で薄く塗る
    mask = grid_z <= threshold
    masked = np.where(mask, 1, np.nan)  # マスクを1にする
    ax.contourf(grid_x, grid_y, masked, levels=[0.5, 1.5], colors=["#a8e6a3"], alpha=0.3)

    # 装飾
    ax.set_title(f"KNN-based AD (k={k}, threshold={threshold})")
    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.legend()
    st.pyplot(fig7)


def shap_display():
    # Streamlit UI
    st.title("🍷 Wine Data - Random Forest 回帰")

    standardize = st.checkbox("標準化する", value=True)

    # データ取得
    X, y = get_wine_data(task="regression", standardize=standardize)

    # データ分割
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # モデル構築と学習
    model = RandomForestRegressor(random_state=42)
    model.fit(X_train, y_train)

    # 予測
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # 評価
    mse_train = mean_squared_error(y_train, y_train_pred)
    r2_train = r2_score(y_train, y_train_pred)
    mse_test = mean_squared_error(y_test, y_test_pred)
    r2_test = r2_score(y_test, y_test_pred)

    # 評価指標の表示
    st.subheader("📈 評価指標")
    st.markdown(f"**訓練データ** - MSE: {mse_train:.3f}, R²: {r2_train:.3f}")
    st.markdown(f"**テストデータ** - MSE: {mse_test:.3f}, R²: {r2_test:.3f}")

    # yyプロット（訓練＋テスト）
    st.subheader("🔍 yyプロット（実測値 vs 予測値）")
    fig, ax = plt.subplots()
    ax.scatter(y_train, y_train_pred, label="train", color="blue", alpha=0.6)
    ax.scatter(y_test, y_test_pred, label="test", color="red", alpha=0.6)
    ax.plot([min(y), max(y)], [min(y), max(y)], 'k--', lw=2)
    ax.set_xlabel("actural")
    ax.set_ylabel("predict")
    ax.set_title("yyplot")
    ax.legend()
    st.pyplot(fig)

    # SHAP解析
    explainer = shap.Explainer(model, X_train)
    shap_values = explainer(X_train, check_additivity=False)

    # データ選択
    st.subheader("解析したいデータポイントを選んでください")
    index = st.slider("サンプル番号", 0, len(X_train) - 1, 0)

    # 特徴量表示
    st.write("### 入力特徴量と値")
    selected_sample = X_train.iloc[index]
    selected_text = "\n".join([f"{col} ({selected_sample[col]:.2f})" for col in X_train.columns])
    st.code(selected_text)

    # 予測と実測
    pred = model.predict([selected_sample])[0]
    true = y_train.iloc[index]
    st.markdown(f"**実測値（Color Intensity）**: {true:.3f}")
    st.markdown(f"**予測値**: {pred:.3f}")

    # force風 summary plot
    st.subheader("特徴量ごとのSHAP値（予測への貢献度）")
    fig, ax = plt.subplots()
    shap.plots.bar(shap_values[index], show=False)
    st.pyplot(fig)

    # SHAP summary plot
    st.subheader("① Summary Plotで全体の傾向を見る")
    fig_summary, ax_summary = plt.subplots()
    shap.plots.beeswarm(shap_values, show=False)
    st.pyplot(fig_summary)

    # SHAP importance bar plot
    st.subheader("② 平均絶対値で特徴量の重要度を見る")
    shap_mean = np.abs(shap_values.values).mean(axis=0)
    shap_importance = pd.DataFrame({
        "Feature": X_train.columns,
        "Mean Absolute SHAP Value": shap_mean
    }).sort_values(by="Mean Absolute SHAP Value", ascending=True)

    fig_bar, ax_bar = plt.subplots()
    shap_importance.plot.barh(x="Feature", y="Mean Absolute SHAP Value", color="red", ax=ax_bar, legend=False)
    ax_bar.set_title("Feature Importance based on Mean Absolute SHAP Values")
    st.pyplot(fig_bar)

    # ▼ ユーザー入力で変数選択
    st.subheader("変数を選んで依存関係を表示")
    feature_list = list(X.columns)
    main_feature = st.selectbox("① 横軸にする特徴量（SHAP対象）", feature_list, index=feature_list.index("flavanoids"))
    color_feature = st.selectbox("② 色分けする特徴量（相互作用を見る）", feature_list, index=feature_list.index("alcohol"))

    # プロット生成
    fig, ax = plt.subplots()
    shap.dependence_plot(
        main_feature,
        shap_values.values,
        X_train,
        interaction_index=color_feature,
        ax=ax,
        show=False
    )
    st.pyplot(fig)


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "Interpretability"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()