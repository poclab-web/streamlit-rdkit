import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic import stmolblock
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

from logic.FingerPrint import morgan_fingerprint
from logic import sascore


# アプリの定義

def rdkit_fr_descriptor_display():
    st.title('RDKit Descriptor Fragments Description 😀')
    st.text("データとしてはTCIで購入可能な約30000分子を脱塩処理した後にDescriptorを取得しています")

    with st.expander("Descriptorの内容"):
        st.title('RDKit Descriptor Description 😀')
        df = pd.read_csv("data/descriptors_name.csv", encoding='shift_jis')
        st.dataframe(df, 2000, 4000)

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

def morgan_fingerprint_display():
    st.title('RDKit + MorganFingerPrint 😀')

    smiles = st.text_input('SMILESを入力', 'CC(=O)C')
    radius = st.number_input('radiusを入力', value = 2)
    nBits = st.number_input('Bitsを入力', value = 2048)
    st.text("入力情報を表示")
    st.code(smiles)
    st.code(radius)
    st.code(nBits)

    img, bitI_morgan, mol = morgan_fingerprint(smiles)

    st.text("1となっているbitの合計数")
    st.code(len(bitI_morgan.keys()))

    st.text("1となっている部分構造を表示")
    st.image(img)

    st.text("分子の構造を表示")

    for atom in mol.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

    img2 = Draw.MolToImage(mol)
    st.image(img2)

    st.text("それぞれのbitに入っている情報を表示 IDを表示させた後 （中心原子の番号, radius) を表示させる")
    for bit, value in bitI_morgan.items():
        st.text(bit)
        st.text(value)

def sascore_display():
    st.title('RDKit + sascore 😀')

    search_smiles = st.text_input('SMILESを入力', 'Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]')

    st.text("smilesを表示")
    st.code(search_smiles)

    try:
        mol = Chem.MolFromSmiles(search_smiles)
        st.header("sascore")
        st.text(sascore.calculateScore(mol))
        st.text("合成のしやすさを、1（易）から10（難）まで値で算出")

        blk = stmolblock.makeblock(search_smiles)
        stmolblock.render_mol(blk)

        st.text("Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions")
        st.code("https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8")

    except:
        st.error("構造を正しく記入してください")

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "Descriptors"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()