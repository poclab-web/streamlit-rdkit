import streamlit as st

def display_sidebar():
    """
    サイドバーにGitHub Issuesリンクや参考文献を追加
    """
    st.sidebar.markdown("### サポート・質問")
    st.sidebar.markdown(
        "[こちら](https://github.com/poclab-web/streamlit-rdkit/issues) で質問やバグを報告してください。"
    )
    st.sidebar.markdown("### 解説本")
    st.sidebar.markdown(
        "[ケモインフォマティクス理論（化学記述子編）](https://zenn.dev/poclabweb/books/chemoinfomatics_theory_descriptor) 1から8までのアプリについて解説"    
    )
    st.sidebar.markdown(
        "[ケモインフォマティクス理論（機械学習編）](https://zenn.dev/poclabweb/books/chemoinfomatics_theory_machinelearning) 9から16までのアプリについて解説"    
    )