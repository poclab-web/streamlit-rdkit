import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, SaltRemover, rdMolDescriptors

import tempfile

class MoleculeAnalyzer:
    def __init__(self, source_type, source):
        """
        Initialize the MoleculeAnalyzer with a source type and source.

        Parameters:
            source_type: str ('sdf' or 'csv')
            source: SDF file object or DataFrame containing SMILES strings
        """
        self.source_type = source_type
        self.source = source
        self.remover = SaltRemover.SaltRemover()

    def analyze(self):
        """
        Analyze the input data (SDF or CSV) and return a DataFrame with relevant details.

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

                # Check for undefined stereo
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                undefined_chiral = any(center[1] == '?' for center in chiral_centers)
                undefined_ez_stereo = any(
                    bond.GetStereo() == Chem.BondStereo.STEREONONE and bond.GetBondType() == Chem.BondType.DOUBLE
                    for bond in mol.GetBonds()
                )
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
