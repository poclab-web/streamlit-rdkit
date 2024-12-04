import streamlit as st

st.set_page_config(page_title="chemoinfomaticstool")

from common.multiapp import MultiApp

from main_pages.GetSmiles import get_smiles
from main_pages.NameSearch import pubchem_search
from main_pages.SmilesSearch import rdkit_smiles_search
from main_pages.MultiSearch import multi_rdkit_smiles_search
from main_pages.DataAnalysis import plotly_analysis
from main_pages.DescriptorsDescription import rdkit_descriptor_description
from main_pages.DescriptorsDescription_fragments import rdkit_fr_descriptor_description
from main_pages.GasteigerCharge import rdkit_charge_analysis
from main_pages.sascore import rdkit_sascore_analysis
from main_pages.FingerPrint import morgan_fingerprint, rdkit_fingerprint
from main_pages.SmartsSearch import rdkit_smarts_search

st.title('Chemoinfomatics Tool')

main = MultiApp()

main.add_app("😀SMILESを出力", get_smiles)
main.add_app("🔍単分子分析(pubchempy)", pubchem_search)
main.add_app("😀SMILES検索(RDKit Descriptor)", rdkit_smiles_search)
main.add_app("⭐️SMARTS検索(TCI Compounds)", rdkit_smarts_search)
main.add_app("📑複数分子分析(RDKit Descriptor)", multi_rdkit_smiles_search)
main.add_app("📝Descriptorの内容(RDKit Descriptor)", rdkit_descriptor_description)
main.add_app("📝Fragmentsの内容(RDKit Descriptor)", rdkit_fr_descriptor_description)
main.add_app("☝️morganfingerprint", morgan_fingerprint)
main.add_app("✌️rdkitfingerprint", rdkit_fingerprint)
main.add_app("🔋ChargeAnalysis", rdkit_charge_analysis)
main.add_app("🧮合成難易度(sascore)", rdkit_sascore_analysis)
main.add_app("📊データ解析(plotly)", plotly_analysis)

main.run()