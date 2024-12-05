import streamlit as st
from streamlit_ketcher import st_ketcher

from common import stmolblock

def get_smiles() -> None:
    """
    描画ツールを使って分子構造を作成し、SMILES文字列と3次元構造を表示します。

    ユーザーは構造を描画した後、「Apply」をクリックすることで、
    SMILES文字列と3次元構造を確認できます。

    Raises:
        Exception: SMILESが入力されなかった場合にエラーメッセージを表示します。
    """
    st.write("構造を書いてApplyをクリックするとSMILESと３次元構造が表示されます")

    # SMILES文字列を取得
    smiles: str = st_ketcher()  # 分子構造エディタからSMILESを取得
    st.write("SMILES")
    st.code(smiles)

    try:
        # SMILES文字列から分子構造のモルブロックを生成
        blk: str = stmolblock.makeblock(smiles)
        stmolblock.render_mol(blk)  # 分子構造を3Dビューでレンダリング
    except Exception:
        # 構造が入力されていない場合のエラーメッセージ
        st.error("構造を入力してください")
