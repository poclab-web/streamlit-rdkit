import inspect
import streamlit as st


def display_code(function, title):
    """
    指定された関数のソースコードを折りたたみ表示で表示する。

    Parameters:
        title (str): セクションのタイトル。
        function (function): 表示する対象の関数。
    """
    st.write(f"📜 {title} ")
    try:
        code = inspect.getsource(function)
        st.code(code, language="python")
    except Exception as e:
        st.error(f"コードを取得できませんでした: {e}")
