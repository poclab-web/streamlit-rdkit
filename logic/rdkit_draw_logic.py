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

def smiles_to_data(smiles_list, fields=None):
    """
    Calculates specified properties (structure, molecular weight, molLogP, etc.) from a list of SMILES strings.
    fields: List of output fields (e.g., ["Structure", "MolecularWeight", "molLogP"])
    """
    if fields is None:
        fields = ["Structure", "MolecularWeight", "molLogP"]  # Default: all fields

    data = []
    for smiles in smiles_list:
        entry = {"SMILES": smiles}
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if "Structure" in fields:
                    img = Draw.MolToImage(mol, size=(150, 150))
                    entry["Structure"] = img
                if "MolecularWeight" in fields:
                    entry["MolecularWeight"] = Descriptors.MolWt(mol)
                if "molLogP" in fields:
                    entry["molLogP"] = Descriptors.MolLogP(mol)
            else:
                for f in fields:
                    entry[f] = "Invalid SMILES" if f == "Structure" else "-"
        except Exception as e:
            for f in fields:
                entry[f] = str(e) if f == "Structure" else "-"
        data.append(entry)
    return data
