import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

def create_yyplot(y_true, y_pred):
    """
    実際の値と予測値の散布図 (y-y プロット) を作成し、MatplotlibのFigureを返します。

    Parameters:
        y_true (array-like): 実際の値
        y_pred (array-like): 予測値

    Returns:
        matplotlib.figure.Figure: 作成されたプロット
    """
    fig, ax = plt.subplots()
    ax.scatter(y_true, y_pred, alpha=0.7, edgecolors="k")
    ax.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], "r--", lw=2)
    ax.set_xlabel("Actual")
    ax.set_ylabel("Predicted")
    ax.set_title("yyplot")
    return fig

def calculate_regression_metrics(y_true, y_pred):
    """
    回帰評価指標を計算する関数。

    Parameters:
        y_true (array-like): 実際の値。
        y_pred (array-like): 予測値。

    Returns:
        dict: 計算された指標を含む辞書。
    """
    mae = mean_absolute_error(y_true, y_pred)
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_true, y_pred)
    return {
        "MAE": mae,
        "MSE": mse,
        "RMSE": rmse,
        "R2": r2
    }

