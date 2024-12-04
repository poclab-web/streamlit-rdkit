import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors

def morgan_fingerprint():
    st.title('RDKit + MorganFingerPrint 😀')

    search_smiles = st.text_input('SMILESを入力', 'CC(=O)C')
    radius = st.number_input('radiusを入力', value = 2)
    nBits = st.number_input('Bitsを入力', value = 2048)
    st.text("入力情報を表示")
    st.code(search_smiles)
    st.code(radius)
    st.code(nBits)

    m = Chem.MolFromSmiles(search_smiles)

    bitI_morgan = {}
    AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits, bitInfo=bitI_morgan)
    morgan_turples = ((m, bit, bitI_morgan) for bit in list(bitI_morgan.keys()))
    img = Draw.DrawMorganBits(morgan_turples, molsPerRow=4, legends=['bit: ' + str(x) for x in list(bitI_morgan.keys())])
    st.text("1となっているbitの合計数")
    st.code(len(bitI_morgan.keys()))
    st.text("1となっている部分構造を表示")
    st.image(img)

    st.text("分子の構造を表示")
    for atom in m.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))
    img2 = Draw.MolToImage(m)
    st.image(img2)

    st.text("それぞれのbitに入っている情報を表示　IDを表示させた後 （中心原子の番号, radius) を表示させる")
    for bit, value in bitI_morgan.items():
        st.text(bit)
        st.text(value)

def rdkit_fingerprint():
    st.title('RDKitFingerPrint 😀')

    search_smiles = st.text_input('SMILESを入力', 'Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]')
    m = Chem.MolFromSmiles(search_smiles)
    img = Draw.MolToImage(m)
    st.image(img)
    maxPath = st.number_input("maxPathを入力", value=5)
    fpSize = st.number_input("fpSizeを入力", value=2048)

    # fingerprintの取得
    bitI_rdkit = {}
    fp_rdkit = Chem.RDKFingerprint(m, maxPath=maxPath, fpSize=fpSize, bitInfo=bitI_rdkit)
    bit_number = st.select_slider('表示させたいbitを指定', options= list(bitI_rdkit))
    img2 = Draw.DrawRDKitBit(m, bit_number, bitI_rdkit)
    st.image(img2)
    display_fragments_number = st.slider('一度に表示させたいフラグメントの数を指定(小さい番号順に表示)', 1, len(bitI_rdkit), 12)
    rdkit_turples = ((m, bit, bitI_rdkit) for bit in list(bitI_rdkit.keys())[:display_fragments_number])
    img3 = Draw.DrawRDKitBits(rdkit_turples, molsPerRow=4, legends=['bit: ' + str(x) for x in list(bitI_rdkit.keys())[:display_fragments_number]])
    st.image(img3)

