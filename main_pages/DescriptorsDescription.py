import streamlit as st
import pandas as pd

# @st.cache_data
# def convert_df(df):
#    return df.to_csv().encode('utf-8')

def rdkit_descriptor_description():
    st.title('RDKit Descriptor Description 😀')
    df = pd.read_csv("data/descriptors_name.csv", encoding='shift_jis')
    st.dataframe(df, 2000, 4000)






