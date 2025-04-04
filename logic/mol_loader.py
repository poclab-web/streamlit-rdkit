import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

class MoleculeDataLoader:
    def __init__(self):
        pass

    @staticmethod
    def load_smiles_file(file_path, delimiter=" ", has_name=True):
        """
        Load SMILES data from a plain text file.

        Args:
            file_path (str): Path to the SMILES file.
            delimiter (str): Delimiter used in the file (default: space).
            has_name (bool): Whether the file includes a name column (default: True).

        Returns:
            pd.DataFrame: DataFrame with SMILES strings and RDKit Molecule objects.
        """
        data = []
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    parts = line.strip().split(delimiter)
                    if has_name and len(parts) >= 2:  # SMILES + name
                        smiles, name = parts[0], parts[1]
                    elif not has_name:  # SMILES only
                        smiles, name = parts[0], None
                    else:
                        continue

                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        entry = {"SMILES": smiles, "Mol": mol}
                        if name:
                            entry["Name"] = name
                        data.append(entry)
        except Exception as e:
            print(f"Error reading SMILES file: {e}")
        return pd.DataFrame(data)

    @staticmethod
    def load_csv(file_path, separator='\t', smiles_col='SMILES'):
        """
        Load dataset from a CSV file and add RDKit Molecule column.

        Args:
            file_path (str): Path to the CSV file.
            separator (str): Delimiter used in the file (default: tab).
            smiles_col (str): Name of the column containing SMILES strings.

        Returns:
            pd.DataFrame: DataFrame with SMILES strings and RDKit Molecule objects.
        """
        try:
            data = pd.read_csv(file_path, sep=separator)
            PandasTools.AddMoleculeColumnToFrame(data, smilesCol=smiles_col, molCol='Molecule', includeFingerprints=False)
            return data
        except Exception as e:
            print(f"Error loading CSV file: {e}")
            return None

    @staticmethod
    def load_excel(file_path, smiles_column='SMILES', sheet_name='Sheet1'):
        """
        Load dataset from an Excel file and add RDKit Molecule column.

        Args:
            file_path (str): Path to the Excel file.
            smiles_column (str): Column name containing SMILES strings.
            sheet_name (str): Sheet name to load data from.

        Returns:
            pd.DataFrame: DataFrame with SMILES strings and RDKit Molecule objects.
        """
        try:
            data = pd.read_excel(file_path, sheet_name=sheet_name)
            PandasTools.AddMoleculeColumnToFrame(data, smiles_column=smiles_column, molCol='Molecule')
            return data
        except Exception as e:
            print(f"Error loading Excel file: {e}")
            return None
