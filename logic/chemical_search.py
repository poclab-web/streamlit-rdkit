import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from molvs import Standardizer
import pandas as pd

# CSVファイルの読み込み
def load_tci_data(file_path):
    try:
        data = pd.read_csv(file_path)
        return data
    except Exception as e:
        st.error(f"データの読み込み中にエラーが発生しました: {e}")
        return None

data_file = 'data/TCI_smiles.csv'
tci_data = load_tci_data(data_file)

if tci_data is not None:
    smiles_list = tci_data['smiles'].tolist()
else:
    st.stop()

# 標準化関数（塩の除去）
def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    standardizer = Standardizer()
    standardized_mol = standardizer.standardize(mol)
    return Chem.MolToSmiles(standardized_mol)

# 完全一致検索関数
def search_exact_match(query_smiles, smiles_list, ignore_stereo=False, include_salts=False):
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        return []

    if not include_salts:
        # クエリを標準化（塩除去）
        query_smiles = standardize_smiles(query_smiles)

    matched_smiles = []
    for smiles in smiles_list:
        target_smiles = smiles

        # 標準化（塩除去）
        if not include_salts:
            target_smiles = standardize_smiles(smiles)

        query_mol = Chem.MolFromSmiles(query_smiles)
        target_mol = Chem.MolFromSmiles(target_smiles)

        if not target_mol:
            continue

        # 立体情報を無視
        if ignore_stereo:
            query_smiles = Chem.MolToSmiles(query_mol, isomericSmiles=False)
            target_smiles = Chem.MolToSmiles(target_mol, isomericSmiles=False)

        if query_smiles == target_smiles:
            matched_smiles.append(smiles)

    return matched_smiles