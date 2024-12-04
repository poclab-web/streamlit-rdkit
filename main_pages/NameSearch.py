import pandas as pd
import streamlit as st

import pubchempy as pcp


def pubchem_search():
    st.title('Pubchempy + Py3DMOL 😀')

    compound_name = st.text_input('名前を入力(cas番号, 慣用名, IUPAC名など)', 'acetone')

    try:
        # PubChemからデータを取得
        df = pcp.get_compounds(compound_name, 'name', as_dataframe=True)

        # データが存在するかを確認
        if df.empty:
            raise ValueError("該当するデータがありません。")

        # canonical_smiles を表示
        st.markdown("### canonical_smiles")
        st.code(df.at[df.index.values[0], "canonical_smiles"])

        # CID を表示
        st.markdown("### CID number")
        st.code(df.index.values[0])

        # 各列のデータを表示
        for index, data in df.items():
            st.markdown("### " + index)  # 列名を表示
            st.code("\n".join(map(str, data.tolist())))  # 列の全データを表示

    except ValueError as ve:
        st.error(str(ve))  # データがない場合のエラーメッセージ
    except Exception as e:
        st.error(f"エラーが発生しました: {str(e)}")  # その他のエラー



