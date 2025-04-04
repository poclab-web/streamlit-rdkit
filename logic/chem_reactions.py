from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

class ReactionSmilesParser:
    def __init__(self, reaction_smiles):
        """
        REACTION SMILESを受け取り、クラス内で処理する準備をする。
        
        :param reaction_smiles: リストまたは単一のREACTION SMILES文字列
        """
        if isinstance(reaction_smiles, str):
            self.reaction_smiles_list = [reaction_smiles]
        elif isinstance(reaction_smiles, list):
            self.reaction_smiles_list = reaction_smiles
        else:
            raise ValueError("Input must be a string or a list of strings.")

    def parse_reaction_smiles(self, reaction_smiles):
        """
        単一のREACTION SMILESをReactants, Reagents, Productsに分解。

        :param reaction_smiles: REACTION SMILES文字列
        :return: Reactants, Reagents, Productsのタプル
        """
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            raise ValueError("Invalid reaction SMILES format. Expected 'Reactants>Reagents>Products'.")

        reactants = parts[0].split('.')
        reagents = parts[1].split('.') if parts[1] else []
        products = parts[2].split('.')

        return reactants, reagents, products

    def to_dataframe(self):
        """
        REACTION SMILESリストをReactants, Reagents, Products列を持つデータフレームに変換。

        :return: pandas.DataFrame
        """
        data = []
        for smiles in self.reaction_smiles_list:
            try:
                reactants, reagents, products = self.parse_reaction_smiles(smiles)
                data.append({
                    "Reactants": reactants,
                    "Reagents": reagents,
                    "Products": products
                })
            except ValueError as e:
                print(f"Error parsing SMILES: {smiles}. {e}")

        return pd.DataFrame(data)
    
    def draw_molecules(self, smiles_list):
        """
        SMILESリストを描画する。

        :param smiles_list: 分子のSMILES文字列のリスト
        :return: PIL Imageオブジェクト
        """
        mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles)]
        return Draw.MolsToImage(mols, subImgSize=(200, 200))
