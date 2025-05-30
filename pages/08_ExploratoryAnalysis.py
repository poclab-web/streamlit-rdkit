import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar
import pandas as pd
import plotly.express as px

from logic.mol_loader import MoleculeDataLoader  

@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')

def plotly_analysis_display():
    st.title('Plotly plot 😀')

    # サンプルデータのセクション
    st.header("サンプルデータ")
    if "selected_dataset" not in st.session_state:
        st.session_state["selected_dataset"] = "solubility"
    if "sample_data" not in st.session_state:
        st.session_state["sample_data"] = None

    dataset_name = st.selectbox(
        "データセットを選択してください",
        ["solubility", "NMR", "food", "qssr"],
        key="selected_dataset"
    )

    if st.button("Load Sample Data"):
        try:
            if dataset_name == "solubility":
                file_path = 'data/curated-solubility-dataset.tab'
                st.session_state["sample_data"] = pd.read_csv(file_path, sep='\t')
            elif dataset_name == "NMR":
                file_path = 'data/NMRshiftDB2_CHOonly_no_missing.xlsx'
                st.session_state["sample_data"] = pd.read_excel(file_path)
            elif dataset_name == "food":
                file_paths = ['data/TasteDB.smi', 'data/FragranceDB.smi']
                data_frames = [MoleculeDataLoader.load_smiles_file(file) for file in file_paths]
                st.session_state["sample_data"] = pd.concat(data_frames, ignore_index=True)
            elif dataset_name == "qssr":
                st.warning("QSSRデータセットのロード機能は未実装です。")
                st.session_state["sample_data"] = None
            else:
                st.error("不明なデータセットが選択されました。")
                st.session_state["sample_data"] = None
            if st.session_state["sample_data"] is not None:
                st.write(f"### サンプルデータ: {dataset_name}")
                st.dataframe(st.session_state["sample_data"])
        except Exception as e:
            st.error(f"サンプルデータのロード中にエラーが発生しました: {e}")

    uploaded_file = st.file_uploader("csvファイルをアップロードしてください")
    test_df = pd.read_csv("data/soac.csv")
    test = convert_df(test_df)

    st.download_button(
        "example csvのDownload",
        test,
        "example.csv",
        "text/csv",
        key='download-csv'
    )

    # 解析対象データの選択
    df = None
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.write("アップロードデータ")
        st.dataframe(df)
    elif st.session_state["sample_data"] is not None:
        df = st.session_state["sample_data"]

    if df is not None:
        X = st.selectbox("select X", df.columns.values.tolist())
        Y = st.selectbox("select Y", df.columns.values.tolist())
        fig = px.scatter(df, x=X, y=Y)
        st.plotly_chart(fig, use_container_width=True)

def filter_display():
    st.header("相関係数の表示")

    # サンプルデータ選択UIを追加
    if "selected_dataset" not in st.session_state:
        st.session_state["selected_dataset"] = "solubility"
    if "sample_data" not in st.session_state:
        st.session_state["sample_data"] = None

    dataset_name = st.selectbox(
        "サンプルデータセットを選択してください（相関係数表示用）",
        ["solubility", "NMR", "food", "qssr"],
        key="corr_selected_dataset"
    )

    load_sample = st.button("サンプルデータをロード（相関係数用）", key="corr_load_sample")
    if load_sample:
        try:
            if dataset_name == "solubility":
                file_path = 'data/curated-solubility-dataset.tab'
                st.session_state["sample_data"] = pd.read_csv(file_path, sep='\t')
            elif dataset_name == "NMR":
                file_path = 'data/NMRshiftDB2_CHOonly_no_missing.xlsx'
                st.session_state["sample_data"] = pd.read_excel(file_path)
            elif dataset_name == "food":
                file_paths = ['data/TasteDB.smi', 'data/FragranceDB.smi']
                data_frames = [MoleculeDataLoader.load_smiles_file(file) for file in file_paths]
                st.session_state["sample_data"] = pd.concat(data_frames, ignore_index=True)
            elif dataset_name == "qssr":
                st.warning("QSSRデータセットのロード機能は未実装です。")
                st.session_state["sample_data"] = None
            else:
                st.error("不明なデータセットが選択されました。")
                st.session_state["sample_data"] = None
            if st.session_state["sample_data"] is not None:
                st.write(f"### サンプルデータ: {dataset_name}")
                st.dataframe(st.session_state["sample_data"])
        except Exception as e:
            st.error(f"サンプルデータのロード中にエラーが発生しました: {e}")

    df = None
    # サンプルデータまたはアップロードデータを利用
    if "sample_data" in st.session_state and st.session_state["sample_data"] is not None:
        df = st.session_state["sample_data"]
    uploaded_file = st.file_uploader("相関係数を計算するcsvファイルをアップロードしてください", key="corr-upload")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.write("アップロードデータ")
        st.dataframe(df)
    if df is not None:
        numeric_df = df.select_dtypes(include=["number"])
        if numeric_df.empty:
            st.warning("数値データがありません。")
        else:
            st.subheader("数値カラムの相関係数")
            corr = numeric_df.corr()
            st.dataframe(corr)
            fig = px.imshow(corr, text_auto=True, color_continuous_scale="RdBu", title="相関係数ヒートマップ")
            st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("データをロードまたはアップロードしてください。")

if __name__ == "__main__":
    current_category = "ExploratoryAnalysis"
    st.write(f"現在のカテゴリー: {current_category}")
    handle_tabs_for_category(current_category)
    display_sidebar()