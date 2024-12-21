import streamlit as st
import seaborn as sns
import pandas as pd
from sklearn.datasets import load_iris, load_wine


def load_dataset_by_name(name):
    # dataset_name =  ["Titanic", "Iris", "Penguins", "Wine"]
    match name:
        case "Titanic":
            return sns.load_dataset('titanic')
        case "Iris":
            iris = load_iris()
            df = pd.DataFrame(iris.data, columns=iris.feature_names)
            df['species'] = pd.Categorical.from_codes(iris.target, iris.target_names)
            return df
        case "Penguins":
            return sns.load_dataset('penguins')
        case "Wine":
            wine = load_wine()
            df = pd.DataFrame(wine.data, columns=wine.feature_names)
            df['target'] = wine.target
            return df
        case _:
            st.error("データセットが見つかりません！")
            return pd.DataFrame()
