import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

def draw_molecule_2d(smiles):
    """
    SMILESから分子の2D構造を描画
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("無効なSMILESが入力されました")
    
    img = Draw.MolToImage(mol, size=(500, 500))
    return img