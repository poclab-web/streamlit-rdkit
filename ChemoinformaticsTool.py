import streamlit as st
import pandas as pd

from utils.yaml_loader import load_yaml
from utils.sidebar import display_sidebar

# YAMLからアプリ定義をロード
app_definitions = load_yaml("app_definitions.yaml")

# ページ設定
st.set_page_config(page_title="Chemoinformatics Tool", layout="wide")

# メインページタイトル
st.title("Chemoinformatics Tool")

st.warning('現在、このアプリは工事中のものが多いです。2025年夏頃にver1.0の完成を目指して改修中です。', icon="⚠️")

st.write("このサイトは、ケモインフォマティクスの簡単なコードについて学ぶためのサイトです")

# カテゴリー一覧を取得
categories = app_definitions["categories"]

# カテゴリーを表形式で表示
st.markdown("### カテゴリー一覧")
category_data = [
    {
        "カテゴリー": f'<a href="/{category_key}">{category_info["title"]}</a>',
        "説明": category_info["description"]
    }
    for category_key, category_info in categories.items()
]
category_df = pd.DataFrame(category_data)

# 表形式で表示 (HTMLリンクを許可する)
st.write(category_df.to_html(escape=False, index=False), unsafe_allow_html=True)

# サイドバーを表示
display_sidebar()
