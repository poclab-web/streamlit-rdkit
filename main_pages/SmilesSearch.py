import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Descriptors

from logic import stmolblock


def rdkit_smiles_search():
    st.title('RDKit + Py3DMOL 😀')

    search_smiles = st.text_input('SMILESを入力', 'CC(=O)C')

    st.text("smilesを表示")
    st.code(search_smiles)

    try:
        blk=stmolblock.makeblock(search_smiles)
        stmolblock.render_mol(blk)

        df = pd.DataFrame({
        'SMILES' :[search_smiles]
        })

        st.subheader('RDKit Desctiptor')
        for i, j in Descriptors.descList:
            try:
                mol = df["SMILES"].map(Chem.MolFromSmiles)
                df[i] = mol.map(j)
                st.write(i)
                st.write(df[i].values[0])
            except:
                pass
    except:
        st.error("構造を正しく記入してください")

