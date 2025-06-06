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

def calculate_confusion_matrix(y_true, y_pred):
    """
    分類問題の混同行列を計算する関数。

    Parameters:
        y_true (array-like): 実際のクラスラベル。
        y_pred (array-like): 予測されたクラスラベル。

    Returns:
        np.ndarray: 混同行列。
    """
    from sklearn.metrics import confusion_matrix
    return confusion_matrix(y_true, y_pred)

def calculate_classification_metrics(y_true, y_pred):
    """
    分類評価指標を計算する関数。

    Parameters:
        y_true (array-like): 実際のクラスラベル。
        y_pred (array-like): 予測されたクラスラベル。

    Returns:
        dict: 計算された指標を含む辞書。
    """
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average='weighted')
    recall = recall_score(y_true, y_pred, average='weighted')
    f1 = f1_score(y_true, y_pred, average='weighted')
    mcc = matthews_corrcoef(y_true, y_pred)
    
    return {
        "Accuracy": accuracy,
        "Precision": precision,
        "Recall": recall,
        "F1 Score": f1,
        "MCC": mcc
    }