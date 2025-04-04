from rdkit import Chem
from rdkit.Chem import DataStructs, Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdFMCS
import numpy as np

def normalized_levenshtein_distance(s1: str, s2: str) -> float:
    """
    Calculate the normalized Levenshtein distance between two strings.

    Args:
        s1 (str): First string.
        s2 (str): Second string.

    Returns:
        float: Normalized Levenshtein distance.
    """
    if len(s1) < len(s2):
        s1, s2 = s2, s1

    if len(s2) == 0:
        return 1.0

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    lev_distance = previous_row[-1]
    max_len = max(len(s1), len(s2))
    normalized_distance = lev_distance / max_len
    return normalized_distance

def calculate_descriptor_distance(mol1, mol2) -> float:
    """
    Calculate the Euclidean distance between two molecules based on RDKit descriptors.

    Args:
        mol1: RDKit molecule object for the first compound.
        mol2: RDKit molecule object for the second compound.

    Returns:
        float: Euclidean distance between the descriptor vectors.
    """
    descriptors = [Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors, Descriptors.NumHAcceptors]
    vec1 = np.array([desc(mol1) for desc in descriptors])
    vec2 = np.array([desc(mol2) for desc in descriptors])
    distance = np.linalg.norm(vec1 - vec2)
    return distance

def calculate_mcs_similarity(mol1, mol2) -> float:
    """
    Calculate the similarity between two molecules based on the Maximum Common Substructure (MCS).

    Args:
        mol1: RDKit molecule object for the first compound.
        mol2: RDKit molecule object for the second compound.

    Returns:
        float: Similarity score based on the MCS.
    """
    mcs_result = rdFMCS.FindMCS([mol1, mol2])
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    
    if mcs_mol is None:
        return 0.0
    
    mcs_size = mcs_mol.GetNumAtoms()
    mol1_size = mol1.GetNumAtoms()
    mol2_size = mol2.GetNumAtoms()
    
    similarity = (2 * mcs_size) / (mol1_size + mol2_size)
    return similarity

def calculate_similarity(smiles1: str, smiles2: str) -> dict:
    """
    Calculate the similarity between two compounds given their SMILES strings using different methods.

    Args:
        smiles1 (str): SMILES string of the first compound.
        smiles2 (str): SMILES string of the second compound.

    Returns:
        dict: Dictionary containing similarity scores from different methods.
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        raise ValueError("Invalid SMILES string provided.")
    
    # Fingerprint similarity
    fp1 = FingerprintMols.FingerprintMol(mol1)
    fp2 = FingerprintMols.FingerprintMol(mol2)
    fingerprint_similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    
    # Normalized Levenshtein distance
    norm_lev_distance = normalized_levenshtein_distance(smiles1, smiles2)
    
    # Descriptor-based Euclidean distance
    descriptor_distance = calculate_descriptor_distance(mol1, mol2)
    
    # MCS similarity
    mcs_similarity = calculate_mcs_similarity(mol1, mol2)
    
    return {
        "fingerprint_similarity": fingerprint_similarity,
        "normalized_levenshtein_distance": norm_lev_distance,
        "descriptor_distance": descriptor_distance,
        "mcs_similarity": mcs_similarity
    }

def find_similar_compounds(query_smiles: str, smiles_list: list, threshold: float = 0.7) -> list:
    """
    Find compounds similar to the query SMILES from a list of SMILES strings.

    Args:
        query_smiles (str): SMILES string of the query compound.
        smiles_list (list): List of SMILES strings to search.
        threshold (float): Similarity threshold for filtering results.

    Returns:
        list: List of SMILES strings that are similar to the query compound.
    """
    similar_compounds = []
    for smiles in smiles_list:
        similarity_scores = calculate_similarity(query_smiles, smiles)
        if similarity_scores["fingerprint_similarity"] >= threshold:
            similar_compounds.append((smiles, similarity_scores))
    return similar_compounds
