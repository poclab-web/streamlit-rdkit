import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')

def multi_rdkit_smiles_search():
    st.title('RDKit Descriptors 😀')
    uploaded_file = st.file_uploader("smilesという名前がついているcsvファイルをアップロードしてください", type=['csv'])
    test_df = pd.read_csv("data/smiles.csv")
    test = convert_df(test_df)

    st.download_button(
        "example csvのDownload",
        test,
        "smiles_example.csv",
        "text/csv",
        key='download-csv'
    )



    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        for i, j in Descriptors.descList:
            mol = df["smiles"].map(Chem.MolFromSmiles)
            df[i] = mol.map(j)

        st.write(df)

        csv = convert_df(df)

        st.download_button(
            "Press to Download",
            csv,
            "RDKit_descriptors.csv",
            "text/csv",
            key='download-csv'
        )

    else:
        st.text("example csvの出力")
        test_df = pd.read_csv("data/smiles.csv")
        for i, j in Descriptors.descList:
            mol = test_df["smiles"].map(Chem.MolFromSmiles)
            test_df[i] = mol.map(j)

        st.write(test_df)

