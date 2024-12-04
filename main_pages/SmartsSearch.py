import random

import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw

@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')

def rdkit_smarts_search():
    st.title('RDKit + Py3DMOL 😀')

    # smartsを入力
    st.text("TCIで売られている化合物(数万個)の化合物をSMARTSで検索")
    search_smarts = st.text_input('SMARTSを入力', 'c1cc([Oh])ccc1')
    patt = Chem.MolFromSmarts(search_smarts)

    # imgを表示
    st.text("smartsで読み込んだ画像表示")
    img = Draw.MolsToGridImage([patt])
    st.image(img)
    st.code(search_smarts)

    # TCIデータの読み込み
    df = pd.read_csv("data/TCI_smiles.csv", encoding='shift_jis', index_col=0)

    mols = []
    for smi in df["smiles"]:
        try:
            mol = Chem.MolFromSmiles(smi)
            mols.append(mol)
        except:
            pass

    matches = [mol for mol in mols if mol.HasSubstructMatch(patt)]
    st.markdown("### 合致した構造の数" + str(len(matches)))

    molsPerRow = st.text_input('構造例として１行に表示させる個数', '3')
    subImgSize = (300, 200)
    number = st.text_input('構造例として表示させたい分子数', '6')

    random_matches = random.sample(matches, int(number))
    img2 = Draw.MolsToGridImage(random_matches, molsPerRow=int(molsPerRow))
    st.image(img2)

    smi2 = [Chem.MolToSmiles(mol) for mol in random_matches]
    df2 = pd.DataFrame(smi2, columns=["smiles"])
    st.dataframe(df2)

    if st.button('合致した構造全てのcsvファイルの作成'):
        smiles_list = [Chem.MolToSmiles(mol) for mol in matches]
        df_download = pd.DataFrame(smiles_list, columns=["smiles"])
        df_download_csv = convert_df(df_download)

        st.download_button(
            "smiles.csvのDownload",
            df_download_csv,
            "smiles.csv",
            "text/csv",
            key='download-csv'
        )
    else:
        st.write('Please Click Start Download button!')