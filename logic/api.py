import requests
import pandas as pd
import xml.etree.ElementTree as ET
import concurrent.futures

class ChemicalSearchApp:
    BASE_PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    BASE_CHEBI_URL = "https://www.ebi.ac.uk/chebi/ws/rest"
    BASE_NCI_URL = "https://cactus.nci.nih.gov/chemical/structure"
    BASE_CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"
    BASE_ZINC_URL = "https://zinc.docking.org/substances"
    BASE_WIKIPEDIA_URL = "https://en.wikipedia.org/api/rest_v1/page/summary"

    def __init__(self, query, query_type, timeout=3):
        """
        Initialize the ChemicalSearchApp class.

        Args:
        - query (str): The chemical query (e.g., SMILES, InChIKey, Name).
        - query_type (str): The type of query ("SMILES", "InChIKey", "Name").
        - timeout (int, optional): Timeout duration for API requests (default: 3 seconds).
        """
        self.query = query.strip()
        self.query_type = query_type
        self.timeout = timeout  # Set the timeout based on user input or default to 3 seconds

    def is_supported(self, api_name):
        """Check if the query type is supported by a specific API."""
        compatibility = {
            "PubChem": ["SMILES", "InChIKey", "Name"],
            "ChEBI": ["SMILES", "Name"],
            "NCI Resolver": ["SMILES", "InChIKey", "Name"],
            "ChEMBL": ["SMILES", "InChIKey", "Name"],
            "ZINC Database": ["SMILES", "Name"],
            "Wikipedia": ["Name"]
        }
        return self.query_type in compatibility.get(api_name, [])

    def search_pubchem(self):
        if not self.is_supported("PubChem"):
            return pd.DataFrame([{"Error": "PubChem does not support this query type."}])

        if self.query_type == "SMILES":
            url = f"{self.BASE_PUBCHEM_URL}/smiles/{self.query}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey/JSON"
        elif self.query_type == "InChIKey":
            url = f"{self.BASE_PUBCHEM_URL}/inchikey/{self.query}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey/JSON"
        elif self.query_type == "Name":
            url = f"{self.BASE_PUBCHEM_URL}/name/{self.query}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey/JSON"

        return self.fetch_data(url)

    def search_chebi(self):
        if not self.is_supported("ChEBI"):
            return pd.DataFrame([{"Error": "ChEBI does not support this query type."}])

        if self.query_type == "SMILES":
            url = f"{self.BASE_CHEBI_URL}/smiles/{self.query}"
        elif self.query_type == "Name":
            url = f"{self.BASE_CHEBI_URL}/name/{self.query}"

        return self.fetch_data(url)

    def search_nci(self):
        if not self.is_supported("NCI Resolver"):
            return pd.DataFrame([{"Error": "NCI Resolver does not support this query type."}])

        url = f"{self.BASE_NCI_URL}/{self.query}/all"
        return self.fetch_data(url)

    def search_chembl(self):
        if not self.is_supported("ChEMBL"):
            return pd.DataFrame([{"Error": "ChEMBL does not support this query type."}])

        if self.query_type == "SMILES":
            url = f"{self.BASE_CHEMBL_URL}/molecule?molecule_structures__canonical_smiles={self.query}"
        elif self.query_type == "InChIKey":
            url = f"{self.BASE_CHEMBL_URL}/molecule?molecule_structures__standard_inchi_key={self.query}"
        elif self.query_type == "Name":
            url = f"{self.BASE_CHEMBL_URL}/molecule?molecule_synonyms__molecule_synonym={self.query}"

        return self.fetch_data(url)

    def search_zinc(self):
        if not self.is_supported("ZINC Database"):
            return pd.DataFrame([{"Error": "ZINC Database does not support this query type."}])

        if self.query_type == "SMILES":
            url = f"{self.BASE_ZINC_URL}/?smiles={self.query}"
        elif self.query_type == "Name":
            url = f"{self.BASE_ZINC_URL}/?name={self.query}"

        return self.fetch_data(url)

    def search_wikipedia(self):
        if not self.is_supported("Wikipedia"):
            return pd.DataFrame([{"Error": "Wikipedia does not support this query type."}])

        url = f"{self.BASE_WIKIPEDIA_URL}/{self.query.replace(' ', '_')}"
        return self.fetch_data(url)

    def fetch_data(self, url):
        """Fetch data from the given URL with a timeout."""
        try:
            response = requests.get(url, headers={"Accept": "application/json"}, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return pd.DataFrame([data]) if isinstance(data, dict) else pd.DataFrame(data)
        except requests.exceptions.RequestException as e:
            return pd.DataFrame([{"Error": str(e)}])

    def get_all_dataframes(self):
        """Fetch results from all APIs in parallel and return them as a dictionary."""
        apis = {
            "PubChem": self.search_pubchem,
            "ChEBI": self.search_chebi,
            "NCI Resolver": self.search_nci,
            "ChEMBL": self.search_chembl,
            "ZINC Database": self.search_zinc,
            "Wikipedia": self.search_wikipedia,
        }

        results = {}
        # Use ThreadPoolExecutor for parallel API requests
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Submit all API requests in parallel
            future_to_api = {executor.submit(search_function): api_name for api_name, search_function in apis.items()}
            for future in concurrent.futures.as_completed(future_to_api):
                api_name = future_to_api[future]
                try:
                    results[api_name] = future.result()
                except Exception as e:
                    results[api_name] = pd.DataFrame([{"Error": str(e)}])

        return results