import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import py3Dmol
from streamlit_ketcher import st_ketcher
from logic.stmolblock import makeblock, render_mol

from logic.pubchem_logic import fetch_pubchem_data

# 現在のカテゴリー（手動設定）
current_category = "ChemInfo"  # 正しいカテゴリーキーを指定
st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

# アプリの定義
def get_smiles():
    """分子構造を描画し、SMILES形式を出力するアプリ。"""
    smiles = st_ketcher()

    # ユーザーがまだ構造を入力していない場合
    if not smiles:
        st.warning("構造が描画されていません。構造を描画してApplyをクリックしてください。")
        return
    
    # SMILESの表示
    st.write("入力されたSMILES:")
    st.code(smiles)

    try:
        # SMILESから分子構造を生成してレンダリング
        blk = makeblock(smiles)
        render_mol(blk)
        st.code(smiles)
        st.code(blk)

    except Exception as e:
        # ユーザー向けエラーメッセージ
        st.error("3次元構造の描画中にエラーが発生しました。構造が正しいか確認してください。")
        st.error(f"エラー内容: {e}")
        

def pubchem_search():
    """PubChem APIを使った単分子分析アプリ。"""
    compound_name = st.text_input("名前を入力", value="acetone", key="name_input")
    if st.button("検索"):
        try:
            results = fetch_pubchem_data(compound_name)
            st.markdown("### canonical_smiles")
            st.code(results["canonical_smiles"])
            st.markdown("### CID number")
            st.code(results["cid"])
            for index, data in results["data_frame"].items():
                st.markdown(f"#### {index}")
                st.code("\n".join(map(str, data.tolist())))
        except Exception as e:
            st.error(f"エラーが発生しました: {e}")

if __name__ == "__main__":
    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()