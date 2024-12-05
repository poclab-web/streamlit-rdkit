import streamlit as st
import yaml

# YAMLファイルのロード
with open("app_definitions.yaml", "r") as file:
    app_definitions = yaml.safe_load(file)

def display_app_overview():
    """
    アプリケーション一覧ページを動的に生成する関数。
    """
    st.title("アプリケーション一覧")
    st.markdown("以下は、各アプリの説明です。サイドバーで選択して直接アクセスできます。")

    # カテゴリー情報を取得
    categories = app_definitions.get("categories", {})
    apps = app_definitions.get("apps", [])

    # カテゴリーごとに表示
    for category, details in categories.items():
        # カテゴリーセクションのデザイン
        st.markdown(f"### **{category}**")
        st.markdown(f"_{details['description']}_")

        # カテゴリーに属するアプリを表示
        category_apps = [app for app in apps if app["category"] == category]

        for app in category_apps:
            # アプリのタイトルと説明を表示
            st.markdown(f"#### {app['title']}")
            st.markdown(f"_{app['description']}_")

        # 区切り線を追加
        st.markdown("---")
