import streamlit as st
from utils.yaml_loader import load_yaml
from utils.app_runner import load_function
from utils.display_code import display_code
import inspect

def handle_tabs_for_category(current_category):
    """
    指定されたカテゴリーに属するアプリをタブ形式で表示し、実行する。

    Parameters:
        current_category (str): 現在のカテゴリーキー。
    """
    # YAMLからアプリ定義をロード
    app_definitions = load_yaml("app_definitions.yaml")

    # カテゴリー情報を取得
    category_info = app_definitions["categories"].get(current_category)
    if not category_info:
        st.error(f"カテゴリー '{current_category}' が見つかりません。")
        st.stop()

    # メインエリアにカテゴリー情報を表示
    st.title(category_info.get("title", ""))
    st.write(category_info.get("description", ""))

    # 現在のカテゴリーに属するアプリを取得

    apps = [
        app for app in app_definitions["apps"]
        if app.get("category") == current_category
    ]

    with st.expander(f"📜 {current_category}の内容一覧"):
        st.write(apps)

    if not apps:
        st.warning("このカテゴリーにはアプリが登録されていません。")
        st.stop()

    # タブを作成（タブ名をユニークにするためインデックスを追加）
    tabs = st.tabs([f"{i+1}: {app['title']}" for i, app in enumerate(apps)])


    # 各タブにアプリを表示
    for i, (tab, app) in enumerate(zip(tabs, apps)):
        with tab:
            # タイトルと説明を表示
            st.subheader(app["title"])
            st.write(app["description"])
      
            with st.expander("以下のアプリのコードの詳細"):
                # コード表示設定（show_codeがTrueの場合）
                if app.get("show_code", False):
                    # Streamlit関数のコード表示
                    if "function" in app:
                        streamlit_function = load_function(app["function"])
                        display_code(streamlit_function, title="Streamlitに表示するためのコード")
                    
                    # logic_functionsを処理
                    if "logic_functions" in app:
                        logic_function = load_function(app["logic_functions"])
                        display_code(logic_function, title="関数のコード")

            # アプリを実行
            try:
                if "function" in app:
                    streamlit_function = load_function(app["function"])
                    streamlit_function()  # 実行
                else:
                    st.warning("このアプリは現在工事中です。")
            except Exception as e:
                st.error(f"アプリの実行中にエラーが発生しました: {e}")



