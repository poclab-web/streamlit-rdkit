import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, SaltRemover, rdMolDescriptors

import tempfile

class MoleculeAnalyzer:
    def __init__(self, source_type, source):
        """
        Initialize the MoleculeAnalyzer with a source type and source.

        Parameters:
            source_type: str ('sdf', 'csv', 'smiles', or 'smiles_list')
            source: SDF file object, DataFrame containing SMILES strings, a single SMILES string, or a list of SMILES strings
        """
        self.source_type = source_type
        self.source = source
        self.remover = SaltRemover.SaltRemover()

    def analyze(self):
        """
        Analyze the input data (SDF, CSV, SMILES, or SMILES list) and return a DataFrame with relevant details.

        Returns:
            pd.DataFrame: Processed DataFrame
            dict: Summary of the analysis
        """
        data = []

        if self.source_type == 'sdf':
            # Save the uploaded SDF file to a temporary file
            with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                tmp_file.write(self.source.read())
                tmp_file_path = tmp_file.name

            supplier = Chem.SDMolSupplier(tmp_file_path)
            for mol in supplier:
                data.extend(self._process_molecule(mol))

        elif self.source_type == 'csv':
            for _, row in self.source.iterrows():
                smiles = row.get("SMILES")
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                data.extend(self._process_molecule(mol, smiles))

        elif self.source_type == 'smiles':
            # Process a single SMILES string
            mol = Chem.MolFromSmiles(self.source)
            data.extend(self._process_molecule(mol, self.source))

        elif self.source_type == 'smiles_list':
            # Process a list of SMILES strings
            for smiles in self.source:
                mol = Chem.MolFromSmiles(smiles)
                data.extend(self._process_molecule(mol, smiles))

        else:
            raise ValueError(f"Unsupported source_type: {self.source_type}")

        df = pd.DataFrame(data)

        # Add duplicate information
        df['IsDuplicate'] = df['SMILES'].duplicated(keep=False)

        # Create summary
        summary = {
            "Total Molecules": len(df),
            "Molecules with Salts": df['HasSalt'].sum(),
            "Molecules with Undefined Stereo": df['HasUndefinedStereo'].sum(),
            "Unique Molecules": len(df['SMILES'].unique()),
            "Duplicate Molecules": df['IsDuplicate'].sum()
        }

        return df, summary

    def analyze_with_normalization(self):
        """
        Analyze the input data with normalization (removing salts and stereochemistry)
        and check for duplicates based on the normalized SMILES.

        Returns:
            pd.DataFrame: Processed DataFrame with normalized SMILES and duplicate information
            dict: Summary of the analysis
        """
        data = []

        if self.source_type == 'sdf':
            with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                tmp_file.write(self.source.read())
                tmp_file_path = tmp_file.name

            supplier = Chem.SDMolSupplier(tmp_file_path)
            for mol in supplier:
                data.extend(self._process_normalized_molecule(mol))

        elif self.source_type == 'csv':
            for _, row in self.source.iterrows():
                smiles = row.get("SMILES")
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                data.extend(self._process_normalized_molecule(mol, smiles))

        elif self.source_type == 'smiles':
            mol = Chem.MolFromSmiles(self.source)
            data.extend(self._process_normalized_molecule(mol, self.source))

        elif self.source_type == 'smiles_list':
            for smiles in self.source:
                mol = Chem.MolFromSmiles(smiles)
                data.extend(self._process_normalized_molecule(mol, smiles))

        else:
            raise ValueError(f"Unsupported source_type: {self.source_type}")

        df = pd.DataFrame(data)

        # Add duplicate information based on normalized SMILES
        df['IsDuplicateNormalized'] = df['NormalizedSMILES'].duplicated(keep=False)

        # Create summary
        summary = {
            "Total Molecules": len(df),
            "Unique Normalized Molecules": len(df['NormalizedSMILES'].unique()),
            "Duplicate Normalized Molecules": df['IsDuplicateNormalized'].sum()
        }

        return df, summary

    def _process_molecule(self, mol, smiles=None):
        """
        Process a single molecule and return relevant details.

        Parameters:
            mol: RDKit Mol object
            smiles: str (optional)

        Returns:
            list of dict: Details of the molecule
        """
        try:
            if mol:
                smiles = smiles or Chem.MolToSmiles(mol)
                stripped_mol = self.remover.StripMol(mol)
                has_salt = mol.GetNumAtoms() != stripped_mol.GetNumAtoms()

                # 不斉点の判定
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                undefined_chiral = any(center[1] == '?' for center in chiral_centers)

                # 二重結合の立体化学情報の判定
                undefined_ez_stereo = any(
                    bond.GetStereo() == Chem.BondStereo.STEREONONE and bond.GetBondType() == Chem.BondType.DOUBLE
                    for bond in mol.GetBonds()
                )

                # 立体化学情報が未定義かどうか
                has_undefined_stereo = undefined_chiral or undefined_ez_stereo

                return [{
                    "SMILES": smiles,
                    "Molecule": mol,
                    "HasSalt": has_salt,
                    "HasUndefinedStereo": has_undefined_stereo
                }]
            else:
                return [{
                    "SMILES": smiles,
                    "Molecule": None,
                    "HasSalt": None,
                    "HasUndefinedStereo": None
                }]
        except Exception:
            return [{
                "SMILES": smiles,
                "Molecule": None,
                "HasSalt": None,
                "HasUndefinedStereo": None
            }]

    def _process_normalized_molecule(self, mol, smiles=None):
        """
        Process a single molecule, normalize it (remove salts and stereochemistry),
        and return relevant details.

        Parameters:
            mol: RDKit Mol object
            smiles: str (optional)

        Returns:
            list of dict: Details of the molecule with normalized SMILES
        """
        try:
            if mol:
                smiles = smiles or Chem.MolToSmiles(mol)
                stripped_mol = self.remover.StripMol(mol)
                Chem.RemoveStereochemistry(stripped_mol)
                normalized_smiles = Chem.MolToSmiles(stripped_mol)

                return [{
                    "SMILES": smiles,
                    "NormalizedSMILES": normalized_smiles,
                    "Molecule": mol
                }]
            else:
                return [{
                    "SMILES": smiles,
                    "NormalizedSMILES": None,
                    "Molecule": None
                }]
        except Exception:
            return [{
                "SMILES": smiles,
                "NormalizedSMILES": None,
                "Molecule": None
            }]

    @staticmethod
    def desalt_smiles(smiles):
        """
        Remove salts and neutralize the molecule from a SMILES string.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None

            # 脱塩処理: 最も大きいフラグメントを取得
            fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            largest_fragment = max(fragments, key=lambda frag: frag.GetNumAtoms(), default=None)
            if not largest_fragment:
                return None

            # 中性化処理: プロトン化または脱プロトン化
            Chem.SanitizeMol(largest_fragment)  # 分子を正規化
            for atom in largest_fragment.GetAtoms():
                # 負電荷を持つ原子にプロトンを追加
                if (atom.GetFormalCharge() < 0):
                    atom.SetFormalCharge(0)
                    atom.UpdatePropertyCache()
                # 正電荷を持つ原子からプロトンを削除
                elif (atom.GetFormalCharge() > 0):
                    atom.SetFormalCharge(0)
                    atom.UpdatePropertyCache()

            # 中性化後のSMILESを返す
            return Chem.MolToSmiles(largest_fragment, isomericSmiles=True)
        except Exception as e:
            return None

    @staticmethod
    def remove_stereo(smiles):
        """
        Remove stereochemistry from a SMILES string.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False)
        except Exception as e:
            return None

def analyze_chirality(smiles):
    """
    分子の不斉点を判定し、不斉点があるのに立体表記がない原子を特定。
    さらに、不斉点における置換基の種類とCIP則における順番を表示。
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, "Invalid SMILES", None
    
    # 不斉点の取得
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    unassigned = [
        (center, label) for center, label in chiral_centers if label == "?"
    ]
    
    # 不斉点の詳細解析
    chiral_details = []
    for center_idx, label in chiral_centers:
        atom = mol.GetAtomWithIdx(center_idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        neighbor_details = []
        
        for nbr_idx in neighbors:
            neighbor_atom = mol.GetAtomWithIdx(nbr_idx)
            neighbor_symbol = neighbor_atom.GetSymbol()
            neighbor_mass = neighbor_atom.GetMass()
            neighbor_details.append((nbr_idx, neighbor_symbol, neighbor_mass))
        
        # CIP順位付けの取得
        if atom.HasProp('_CIPCode'):
            cip_code = atom.GetProp('_CIPCode')
        else:
            cip_code = "Unassigned"
        
        chiral_details.append({
            "center_idx": center_idx,
            "label": label,
            "cip_code": cip_code,
            "substituents": sorted(
                neighbor_details, key=lambda x: x[2], reverse=True
            )  # 質量で暫定的にソート
        })
    
    return mol, chiral_centers, unassigned, chiral_details
