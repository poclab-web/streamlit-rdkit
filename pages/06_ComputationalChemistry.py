import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic.openbabel_utils import smiles_to_3d_with_make3D
from logic.stmolblock import makeblock, render_mol
from rdkit.Chem.Draw import rdMolDraw2D

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
    st.image(img2)

    # 修正箇所: MolDraw2DCairo を使用して直接描画
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)  # 300x300 の画像サイズ
    SimilarityMaps.GetSimilarityMapFromWeights(mol_with_charge, weights=atom_charges, colorMap='bwr', draw2d=drawer)
    drawer.FinishDrawing()

    # 描画結果をバイナリデータとして取得し、Streamlit に表示
    img_data = drawer.GetDrawingText()
    st.image(img_data)

def visualize_smiles_to_3d_with_make3D():
    # Streamlitアプリケーション
    st.title("SMILESから3次元構造を生成・可視化")

    # ユーザー入力
    smiles = st.text_input("SMILESを入力してください", "CCO")  # デフォルトでエタノール

    if st.button("3D構造を生成"):
        with st.spinner("3D構造を生成しています..."):
            # SMILESから3次元構造を生成
            mol_3d = smiles_to_3d_with_make3D(smiles)
            if "Error" in mol_3d or not mol_3d.strip():
                st.error("無効なSMILES形式です。再度入力してください。")
            else:
                render_mol(mol_3d)
                st.code(mol_3d)
    else:
        st.write("SMILES形式を入力し、ボタンを押して分子の3次元構造を生成・可視化してください。")

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ComputationalChemistry"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()