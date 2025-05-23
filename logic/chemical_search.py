import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

from molvs import Standardizer
from molvs.fragment import LargestFragmentChooser

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
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # フラグメントを最大のものだけに限定（塩除去）
        largest = LargestFragmentChooser()
        mol = largest.choose(mol)

        # 明示的に構造の整合性をチェック
        Chem.SanitizeMol(mol)

        # 立体情報などを含んだ標準的なSMILESに変換
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    except Exception as e:
        # 無効な構造はスキップ
        print(f"[standardize_smiles] 無効な構造をスキップ: {smiles} ({e})")
        return None

def search_exact_match(query_smiles, smiles_list, ignore_stereo=False, include_salts=False):
    # クエリの標準化（必要に応じて塩除去）
    if not include_salts:
        query_smiles_std = standardize_smiles(query_smiles)
        if query_smiles_std is None:
            return []
    else:
        query_smiles_std = query_smiles

    try:
        query_mol = Chem.MolFromSmiles(query_smiles_std)
        if query_mol is None:
            return []
        if ignore_stereo:
            Chem.RemoveStereochemistry(query_mol)
    except Exception as e:
        print(f"[クエリ分子読み込み失敗] {query_smiles_std} ({e})")
        return []

    matched_smiles = []
    for smiles in smiles_list:
        try:
            if not include_salts:
                target_smiles_std = standardize_smiles(smiles)
                if target_smiles_std is None:
                    continue
            else:
                target_smiles_std = smiles

            target_mol = Chem.MolFromSmiles(target_smiles_std)
            if target_mol is None:
                continue
            if ignore_stereo:
                Chem.RemoveStereochemistry(target_mol)

            # Molオブジェクトで構造が一致するかをチェック
            if query_mol.HasSubstructMatch(target_mol) and target_mol.HasSubstructMatch(query_mol):
                matched_smiles.append(smiles)

        except Exception as e:
            print(f"[検索中スキップ] SMILES: {smiles} ({e})")
            continue  # 明示的にスキップ

    return matched_smiles