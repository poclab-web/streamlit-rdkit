import streamlit as st
import yaml
from common.multiapp import MultiApp
from common.display_app_overview import display_app_overview

# 各ページの関数をインポート
from main_pages.GetSmiles import get_smiles
from main_pages.NameSearch import pubchem_search
from main_pages.SmilesSearch import rdkit_smiles_search
from main_pages.MultiSearch import multi_rdkit_smiles_search
from main_pages.DataAnalysis import plotly_analysis
from main_pages.DescriptorsDescription import rdkit_descriptor_description
from main_pages.DescriptorsDescription_fragments import rdkit_fr_descriptor_description
from main_pages.GasteigerCharge import rdkit_charge_analysis
from main_pages.sascore import rdkit_sascore_analysis
from main_pages.FingerPrint import morgan_fingerprint, rdkit_fingerprint
from main_pages.SmartsSearch import rdkit_smarts_search

# YAMLファイルのロード
with open("app_definitions.yaml", "r") as file:
    app_definitions = yaml.safe_load(file)

# MultiApp クラスのインスタンスを作成
main = MultiApp()

# 各アプリを登録
for app in app_definitions["apps"]:
    if "function" in app and app["function"] is not None:
        main.add_app(
            title=app["title"],
            func=globals()[app["function"]],
            category=app["category"],
            show_code=app.get("show_code", True)
        )
    else:
        # function がない場合はダミー関数で登録
        def under_construction():
            st.error("このページは現在工事中です。後日ご利用ください。")
        main.add_app(
            title=app["title"],
            func=under_construction,
            category=app["category"],
            show_code=False
        )

# アプリ一覧を登録
main.add_app(
    title="📋 アプリ一覧",
    func=display_app_overview,
    category="アプリ一覧",
    show_code=False
)

# サイドバーにカテゴリー選択を表示
categories = sorted(set(app['category'] for app in main.apps if app['category'] != "アプリ一覧"))
categories.append("アプリ一覧")  # "アプリ一覧" を最後に追加

# カテゴリーの選択
selected_category = st.sidebar.radio(
    "カテゴリーを選択してください",
    options=categories,
    key="category_selector"
)

# カテゴリー内のアプリを表示
selected_app = None
apps_in_category = [app for app in main.apps if app['category'] == selected_category]

selected_title = st.sidebar.radio(
    f"{selected_category} のアプリを選択してください",
    options=[app['title'] for app in apps_in_category],
    key="app_selector"
)

for app in apps_in_category:
    if app['title'] == selected_title:
        selected_app = app


# GitHub Issuesへのリンクを追加
st.sidebar.markdown("---")
st.sidebar.markdown("### サポート・質問")
st.sidebar.markdown(
    "[こちら](https://github.com/poclab-web/streamlit-rdkit/issues) で質問やバグを報告してください。"
)

# メインページに選択したアプリを表示
st.title("Chemoinformatics Tool")

if selected_app:
    # 選択されたアプリの情報を表示
    st.subheader(f"現在選択されているアプリ: {selected_app['title']}")

    # ソースコード表示オプションを配置
    if selected_app.get('show_code', True):
        with st.expander("ソースコードを表示"):
            try:
                import inspect
                code = inspect.getsource(selected_app['function'])
                st.code(code, language="python")
            except Exception as e:
                st.error(f"コードを取得できませんでした: {e}")

    # アプリを実行
    try:
        selected_app['function']()
    except Exception as e:
        st.error(f"アプリの実行中にエラーが発生しました: {e}")

else:
    st.write("左のメニューからカテゴリーとアプリを選択してください。")
