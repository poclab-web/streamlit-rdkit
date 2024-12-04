import streamlit as st
from streamlit_ketcher import st_ketcher

from common import stmolblock

def get_smiles():
    st.write("構造を書いてApplyをクリックするとSMILESと３次元構造が表示されます")
    smiles = st_ketcher()
    st.write("SMILES")
    st.code(smiles)

    try:
        blk = stmolblock.makeblock(smiles)
        stmolblock.render_mol(blk)
    except:
        st.error("構造を入力してください")