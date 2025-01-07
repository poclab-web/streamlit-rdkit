import streamlit as st
from utils.app_runner import load_function
from utils.yaml_loader import load_yaml
from utils.display_code import display_code

current_category = "ChemInfo"  # 正しいカテゴリーキーを指定
st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

function_path = "logic.pubchem_logic.fetch_pubchem_data"

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


for app in apps:
    st.write(app)
    if "logic_functions" in app and app["logic_functions"] is not None:
                    st.write("logic_functionsを動かします。")
                    st.write(app["logic_functions"])
                    logic_function = load_function(app["logic_functions"])
                    display_code(logic_function, title="関数のコード")








# # 関数をロード
# loaded_function = load_function(function_path)



# if loaded_function:
#     # 関数が正常にロードされた場合、利用可能
#     result = loaded_function()
#     st.write(f"結果: {result}")
# else:
#     st.error("関数をロードできませんでした。")
