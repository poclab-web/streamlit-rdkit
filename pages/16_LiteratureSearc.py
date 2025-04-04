import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

# アプリの定義

import streamlit as st
from logic.ChemRxiv_api import fetch_chemrxiv_data, parse_chemrxiv_data  # ロジック部分をインポート


def chemrxiv_search_display():
    # Streamlitのインターフェース
    st.title("ChemRxiv 検索ツール")
    st.write("指定したキーワードで [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/public-dashboard) APIを検索します。")

    # ユーザー入力
    keyword = st.text_input("検索キーワードを入力してください", "machinelearning")

    # 検索ボタン
    if st.button("検索"):
        with st.spinner("検索中..."):
            raw_data = fetch_chemrxiv_data(keyword)
            parsed_data = parse_chemrxiv_data(raw_data)

        if "error" in parsed_data:
            st.error(f"エラー: {parsed_data['error']}")
        else:
            st.success(f"検索結果の総数: {parsed_data['totalCount']}")
            for item in parsed_data["items"]:
                st.subheader(item["title"])
                
                # DOIをリンクとして表示
                if item['doi'] != "DOIなし":
                    doi_url = f"https://doi.org/{item['doi']}"
                    st.write(f"**DOI:** [{item['doi']}]({doi_url})")
                else:
                    st.write("**DOI:** DOIなし")

                st.write(f"**著者:** {item['authors']}")
                st.write(f"**公開日:** {item['publishedDate']}")
                st.write(f"**概要:** {item['abstract']}")
                st.write(f"**キーワード:** {item['keywords']}")
                
                # ライセンスリンクの表示
                if item['license_url'] != "URLなし":
                    st.write(f"**ライセンス:** [{item['license_name']}]({item['license_url']})")
                else:
                    st.write(f"**ライセンス:** {item['license_name']}")
                
                st.markdown("---")


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "LiteratureSearch"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()