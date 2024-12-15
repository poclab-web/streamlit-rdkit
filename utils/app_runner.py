import streamlit as st
import importlib
from utils.display_code import display_code

def load_function(function_path):
    """
    関数のフルパスを指定して関数オブジェクトをロードする。

    Parameters:
        function_path (str): モジュール名と関数名を含むフルパス (例: "module.submodule.function").

    Returns:
        function: 指定された関数オブジェクト。
    """
    try:
        # モジュール名と関数名を分割
        module_name, function_name = function_path.rsplit(".", 1)
        # モジュールをインポート
        module = importlib.import_module(module_name)
        # 関数を取得
        return getattr(module, function_name)
    except ImportError as e:
        st.error(f"モジュール '{module_name}' をインポートできませんでした: {e}")
    except AttributeError as e:
        st.error(f"モジュール '{module_name}' に関数 '{function_name}' が見つかりません: {e}")
    return None
