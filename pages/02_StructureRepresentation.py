import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import time
import concurrent.futures

from logic.rdkit_draw_logic import draw_molecule_2d
from logic.molecularconverter import MoleculeConverter
from logic.api import ChemicalSearchApp

from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdmolops

import pandas as pd
import matplotlib.pyplot as plt


# アプリの関数定義

# smilesから出力される情報を表示

def DisplayChemicalSearch():
    st.title("Chemical Information Search App")
    st.markdown("""
    This application allows you to search for chemical information using:
    - **SMILES**
    - **InChIKey**
    - **Chemical Name**
    """)

    # Add a table to show API compatibility with descriptions
    st.subheader("API Compatibility Table")
    compatibility_data = {
        "API": ["PubChem", "ZINC Database", "ChEMBL", "Wikipedia", "ChEBI", "NCI Resolver"],
        "Description": [
            "Comprehensive database for chemical compounds and properties.",
            "Database for purchasable compounds and screening libraries.",
            "Database for bioactive molecules and drug discovery.",
            "General chemical compounds and their properties.",
            "Focuses on biologically relevant chemical entities.",
            "Chemical structure resolver with multiple input formats."
        ],
        "SMILES": ["✅ Supported", "✅ Supported", "✅ Supported", "❌ Not Supported", "✅ Supported", "✅ Supported"],
        "InChIKey": ["✅ Supported", "❌ Not Supported", "✅ Supported", "❌ Not Supported", "❌ Not Supported", "✅ Supported"],
        "Name": ["✅ Supported", "✅ Supported", "✅ Supported", "✅ Supported", "✅ Supported", "✅ Supported"],
        "Estimated Compounds": [
            110_000_000,   # PubChem
            50_000_000,    # ZINC Database (revised for purchasable compounds)
            2_000_000,     # ChEMBL
            150_000,       # Wikipedia (Category: "Chemical compounds")
            65_000,        # ChEBI
            "unknown"      # NCI Resolver (Unknown)
        ]
    }
    compatibility_df = pd.DataFrame(compatibility_data)
    st.table(compatibility_df)

    # Default values for SMILES, InChIKey, and Name
    default_values = {
        "SMILES": "CCO",  # Ethanol
        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",  # Ethanol InChIKey
        "Name": "ethanol"  # Ethanol common name
    }

    # User input for query type and query
    query_type = st.selectbox("Select the type of input:", ["SMILES", "InChIKey", "Name"])
    query_input = st.text_input(f"Enter a {query_type}:", value=default_values.get(query_type, ""))

    # User input for timeout duration
    timeout = st.slider("Set the timeout duration (seconds):", min_value=1, max_value=30, value=3)

    if st.button("Search"):
        if query_input:
            # Create instance of ChemicalSearchApp
            app = ChemicalSearchApp(query_input, query_type, timeout=timeout)

            # Fetch data from all APIs
            st.subheader("Results")
            with st.spinner("Fetching data from APIs..."):
                # Fetch data in the fixed order
                api_order = ["PubChem", "ZINC Database", "ChEMBL", "Wikipedia", "ChEBI", "NCI Resolver"]
                all_results = app.get_all_dataframes()

                # Display results in the fixed order
                for api_name in api_order:
                    df = all_results.get(api_name, pd.DataFrame([{"Error": "No data available"}]))
                    st.subheader(f"{api_name} Results")
                    if not df.empty and "Error" in df.columns:
                        st.error(df.iloc[0]["Error"])
                    elif not df.empty:
                        st.dataframe(df)
                    else:
                        st.warning(f"No data returned from {api_name}.")
        else:
            st.warning(f"Please enter a valid {query_type}.")

# 2次元情報について表示
def display_adjacency_matrix ():
    # Streamlitアプリのタイトル
    st.title("RDKitで分子の隣接行列を表示")

    # ユーザー入力：SMILES文字列
    smiles = st.text_input("SMILESを入力してください", "C(C=O)")  # デフォルト値はエタノール

    st.write(f"入力されたSMILES: `{smiles}`")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # インデックス番号付きの分子構造を描画
            st.write("分子構造とインデックス:")
            mol_with_indices = Chem.Mol(mol)
            for atom in mol_with_indices.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
            img = Draw.MolToImage(mol_with_indices, size=(400, 400))
            st.image(img, caption="インデックス番号付きの分子構造")

            # 隣接行列の計算
            adj_matrix = rdmolops.GetAdjacencyMatrix(mol)

            # 隣接行列をDataFrameとして表示
            st.write("#### 隣接行列:")
            adj_df = pd.DataFrame(adj_matrix)
            adj_df.index.name = "Index"
            adj_df.columns = [f"Atom {i}" for i in range(len(adj_matrix))]
            st.write(adj_df)


            # 結合数行列をDataFrameとして表示
            st.write("#### 結合数行列:")
            num_atoms = mol.GetNumAtoms()
            bond_order_matrix = [[0] * num_atoms for _ in range(num_atoms)]
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType()
                # 結合の種類を数値に変換
                if bond_type == Chem.rdchem.BondType.SINGLE:
                    bond_order_matrix[i][j] = 1
                    bond_order_matrix[j][i] = 1
                elif bond_type == Chem.rdchem.BondType.DOUBLE:
                    bond_order_matrix[i][j] = 2
                    bond_order_matrix[j][i] = 2
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    bond_order_matrix[i][j] = 3
                    bond_order_matrix[j][i] = 3
                elif bond_type == Chem.rdchem.BondType.AROMATIC:
                    bond_order_matrix[i][j] = 1.5
                    bond_order_matrix[j][i] = 1.5

            bond_order_df = pd.DataFrame(bond_order_matrix)
            bond_order_df.index.name = "Index"
            bond_order_df.columns = [f"Atom {i}" for i in range(num_atoms)]
            st.write(bond_order_df)

            # 距離行列の計算
            distance_matrix = rdmolops.GetDistanceMatrix(mol)

            # 距離行列をDataFrameとして表示
            st.write("### 距離行列（結合距離基準）")
            st.write("分子内の各原子ペア間の最短距離（結合の数）を算出")
            distance_df = pd.DataFrame(distance_matrix)
            distance_df.index.name = "Atom Index"
            distance_df.columns = [f"Atom {i}" for i in range(distance_matrix.shape[1])]
            st.dataframe(distance_df)

            # SDF形式のMolBlockを取得
            sdf_block = Chem.MolToMolBlock(mol)

            lines = sdf_block.split("\n")
            atom_count = int(lines[3][:3].strip())  # 原子数
            bond_count = int(lines[3][3:6].strip())  # 結合数
            
            # 結合情報を抽出
            bond_data = []
            for i in range(atom_count + 4, atom_count + 4 + bond_count):
                bond_info = lines[i].split()
                bond_data.append({
                    "atom_index_1": int(bond_info[0]) - 1,  # 0-based index
                    "atom_index_2": int(bond_info[1]) - 1,
                    "bond_type": int(bond_info[2])
                })

            # 接続表を抽出して表示
            st.write("### 接続表:")
            st.table(bond_data)

            # SDF形式のファイルを出力
            st.write("### SDF形式のMolBlock:")
            st.code(sdf_block, language="plaintext")

        else:
            st.error("無効なSMILESです。正しい形式で入力してください。")
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")

# 3次元情報について表示
def run_molecule_converter_display():
    st.title("分子変換ツール: SMILES, XYZ, SDF相互変換およびZ-matrix生成")

    converter = MoleculeConverter()
    option = st.radio("入力形式を選択", ["SMILES", "XYZ座標", "SDFファイル"])

    if option == "SMILES":
        smiles = st.text_input("SMILESを入力", "CCO")
        if smiles:
            converter.load_from_smiles(smiles)

    elif option == "XYZ座標":
        xyz_input = st.text_area("XYZ座標を入力", "C 0.0 0.0 0.0\nC 1.4 0.0 0.0\nO 2.8 0.0 0.0")
        if xyz_input:
            xyz_list = [line.split() for line in xyz_input.split("\n")]
            xyz_list = [(p[0], float(p[1]), float(p[2]), float(p[3])) for p in xyz_list]
            converter.load_from_xyz(xyz_list)

    elif option == "SDFファイル":
        uploaded_file = st.file_uploader("SDFファイルをアップロード", type="sdf")
        if uploaded_file:
            with open("temp.sdf", "wb") as f:
                f.write(uploaded_file.getbuffer())
            converter.load_from_sdf("temp.sdf")

    # 出力形式
    output_option = st.radio("出力形式を選択", ["Z-matrix", "XYZ形式", "SDF形式"])

    if output_option == "Z-matrix":
        zmatrix = converter.generate_zmatrix()
        df = pd.DataFrame(
            zmatrix,
            columns=["Element", "Bond Atom", "Bond Length", "Angle Atom", "Bond Angle", "Dihedral Atom", "Dihedral Angle"]
        )
        st.write("### Z-matrix")
        st.dataframe(df.fillna("-"))


    elif output_option == "XYZ形式":
        st.write("### XYZ形式")
        xyz_output = converter.export_xyz()
        st.code(xyz_output, language="plaintext")

    elif output_option == "SDF形式":
        output_filename = "output.sdf"

        # SDFファイルをエクスポート
        converter.export_sdf(output_filename)

        # SDFファイルを画面に表示
        st.write("### SDFファイルの内容")
        with open(output_filename, "r") as f:
            sdf_content = f.read()
            st.text_area("SDF形式の内容", sdf_content, height=300)

        # SDFファイルをダウンロードリンクとして提供
        st.write("### SDFファイルをダウンロード")
        with open(output_filename, "rb") as f:
            st.download_button("SDFファイルをダウンロード", f, file_name="output.sdf")



if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "StructureRepresentation"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()