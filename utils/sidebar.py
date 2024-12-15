import streamlit as st

def display_sidebar():
    """
    サイドバーにGitHub Issuesリンクを追加
    """
    st.sidebar.markdown("### サポート・質問")
    st.sidebar.markdown(
        "[こちら](https://github.com/poclab-web/streamlit-rdkit/issues) で質問やバグを報告してください。"
    )
