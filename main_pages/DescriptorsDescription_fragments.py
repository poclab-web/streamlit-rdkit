import streamlit as st
import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem

# @st.cache_data
# def convert_df(df):
# return df.to_csv().encode('utf-8')

def rdkit_fr_descriptor_description():
    st.title('RDKit Descriptor Fragments Description 😀')
    st.text("データとしてはTCIで購入可能な約30000分子を脱塩処理した後にDescriptorを取得しています")
    df = pd.read_csv("data/TCI_fr.csv", encoding='shift_jis')
    fr_list = df.columns[3:]
    option = st.selectbox('fragmentを選択', fr_list)
    molsPerRow = st.text_input('１行に表示させる個数', '2')
    subImgSize = (300, 200)
    number = st.text_input('表示させたい分子数', '6')
    df2 = df[(df[option] >= 1)]

    mols = []
    smiles_list = []

    for smiles in df2["smiles"].sample(int(number)):
        smiles_list.append(smiles)
        mol = Chem.MolFromSmiles(smiles)
        mols.append(mol)

    st.text("選ばれたfragmentを有する分子を表示")
    img = Draw.MolsToGridImage(mols, molsPerRow=int(molsPerRow), legends = smiles_list, subImgSize = (400, 300))
    st.image(img)






