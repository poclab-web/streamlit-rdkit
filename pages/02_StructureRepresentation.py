import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic.rdkit_draw_logic import draw_molecule_2d


# アプリの関数定義
def Draw2DStructure():
    st.markdown("SMILESを入力すると、分子の2D構造が描画されます。")

    # SMILES入力
    smiles = st.text_input("SMILESを入力してください", value="CCO", key="smiles_input") # デフォルトはエタノール

    if st.button("構造を描画"):
        try:
            img = draw_molecule_2d(smiles)
            st.image(img, caption="2D構造", use_container_width=True)
        except Exception as e:
            st.error(f"エラーが発生しました: {e}")

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "StructureRepresentation"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()