import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors

def morgan_fingerprint(smiles, radius=2, nBits=2048):

    mol = Chem.MolFromSmiles(smiles)

    bitI_morgan = {}
    AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits, bitInfo=bitI_morgan)

    morgan_turples = ((mol, bit, bitI_morgan) for bit in list(bitI_morgan.keys()))
    img = Draw.DrawMorganBits(morgan_turples, molsPerRow=4, legends=['bit: ' + str(x) for x in list(bitI_morgan.keys())])

    return img, bitI_morgan, mol



# TODO
# def rdkit_fingerprint(smiles, maxPath=5, fpSize=2048, bit_number):
#     mol = Chem.MolFromSmiles(smiles)
#     mol_img = Draw.MolToImage(mol)
    
#     bitI_rdkit = {}
#     fp_rdkit = Chem.RDKFingerprint(mol, maxPath=maxPath, fpSize=fpSize, bitInfo=bitI_rdkit)
#     frag_img = Draw.DrawRDKitBit(mol, bit_number, bitI_rdkit)

#     return mol_img, frag_img

# def rdkit_fingerprint_display():
#     st.title('RDKitFingerPrint 😀')

#     search_smiles = st.text_input('SMILESを入力', 'Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]')

#     maxPath = st.number_input("maxPathを入力", value=5)
#     fpSize = st.number_input("fpSizeを入力", value=2048)
#     bit_number = st.select_slider('表示させたいbitを指定', options= list(bitI_rdkit))

#     st.image(mol_img)
#     st.image(frag_img)

#     display_fragments_number = st.slider('一度に表示させたいフラグメントの数を指定(小さい番号順に表示)', 1, len(bitI_rdkit), 12)
#     rdkit_turples = ((mol, bit, bitI_rdkit) for bit in list(bitI_rdkit.keys())[:display_fragments_number])
#     img3 = Draw.DrawRDKitBits(rdkit_turples, molsPerRow=4, legends=['bit: ' + str(x) for x in list(bitI_rdkit.keys())[:display_fragments_number]])
#     st.image(img3)

