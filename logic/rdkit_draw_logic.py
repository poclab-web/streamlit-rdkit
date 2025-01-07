import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

import pandas as pd

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

def smiles_to_data(smiles_list):
    """
    SMILES文字列のリストから構造、分子量、molLogPを計算
    """
    data = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # 分子量とmolLogPを計算
                mol_weight = Descriptors.MolWt(mol)
                mol_logp = Descriptors.MolLogP(mol)
                # 構造の描画
                img = Draw.MolToImage(mol, size=(150, 150))
                data.append({"SMILES": smiles, "構造": img, "分子量": mol_weight, "molLogP": mol_logp})
            else:
                data.append({"SMILES": smiles, "構造": "無効なSMILES", "分子量": "-", "molLogP": "-"})
        except Exception as e:
            data.append({"SMILES": smiles, "構造": str(e), "分子量": "-", "molLogP": "-"})
    return data
