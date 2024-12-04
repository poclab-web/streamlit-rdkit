import streamlit as st

from rdkit import Chem
from common import sascore
from common import stmolblock

def rdkit_sascore_analysis():
    st.title('RDKit + sascore 😀')

    search_smiles = st.text_input('SMILESを入力', 'Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]')

    st.text("smilesを表示")
    st.code(search_smiles)

    try:
        mol = Chem.MolFromSmiles(search_smiles)
        st.header("sascore")
        st.text(sascore.calculateScore(mol))
        st.text("合成のしやすさを、1（易）から10（難）まで値で算出")
        blk=stmolblock.makeblock(search_smiles)
        stmolblock.render_mol(blk)

        st.text("Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions")
        st.code("https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8")

    except:
        st.error("構造を正しく記入してください")


