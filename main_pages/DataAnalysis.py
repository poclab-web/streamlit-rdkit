import streamlit as st
import pandas as pd
import plotly.express as px


@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')

def plotly_analysis():
    st.title('Plotly plot 😀')
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

    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)

        st.dataframe(df)

        X = st.selectbox("select X", df.columns.values.tolist())
        Y = st.selectbox("select Y", df.columns.values.tolist())

        fig = px.scatter(df, x= X, y = Y)
        st.plotly_chart(fig, use_container_width=True)







