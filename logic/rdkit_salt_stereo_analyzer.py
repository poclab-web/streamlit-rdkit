import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, SaltRemover
import streamlit as st

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
