import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar



# アプリの定義

import streamlit as st

from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps

def gasteiger_charge_desplay():
    st.title('RDKit + GasteigerCharge 😀')

    search_smiles = st.text_input('SMILESを入力', 'c1ncncc1C(=O)[O-]')
    mol = Chem.MolFromSmiles(search_smiles)
    mol = Chem.AddHs(mol)
    mol_with_charge = Chem.Mol(mol)

    AllChem.ComputeGasteigerCharges(mol_with_charge)
    atom_charges = [float(mol_with_charge.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(mol_with_charge.GetNumAtoms())]

    for atom in mol_with_charge.GetAtoms():
        lbl = str(round(atom.GetDoubleProp("_GasteigerCharge"), 3))
        atom.SetProp('atomNote', lbl)

    img2 = Draw.MolToImage(mol_with_charge)
    st.image(img2, use_container_width=True)

    fig = SimilarityMaps.GetSimilarityMapFromWeights(mol_with_charge, weights=atom_charges, colorMap='bwr', size=(300, 300))
    fig.savefig('sample.png', bbox_inches='tight')
    image = Image.open('sample.png')
    st.image(image, use_container_width=True)


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ComputationalChemistry"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()