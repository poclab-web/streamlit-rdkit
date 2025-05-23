import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import py3Dmol
from streamlit_ketcher import st_ketcher
import pandas as pd
from rdkit import Chem

from logic.stmolblock import makeblock, render_mol
from logic.pubchem_logic import fetch_pubchem_data
from logic.rdkit_draw_logic import smiles_to_data, draw_molecule_2d

# 現在のカテゴリー（手動設定）
current_category = "ChemInfo"  # 正しいカテゴリーキーを指定
st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

# アプリの定義
def get_smiles():
    """分子構造を描画し、SMILES形式を出力するアプリ。"""
    smiles = st_ketcher()

    # ユーザーがまだ構造を入力していない場合
    if not smiles:
        st.warning("構造が描画されていません。構造を描画してApplyをクリックしてください。")
        return
    
    # SMILESの表示
    st.write("入力されたSMILES:")
    st.code(smiles)

    try: 
        mol = Chem.MolFromSmiles(smiles)

        st.write("InChi")
        st.code(Chem.MolToInchi(mol))

        st.write("InChiKey")
        st.code(Chem.MolToInchiKey(mol))  # InChIKeyを表示

        col1, col2 = st.columns(2)

        # Display 2D structure in the first column
        with col1:
            st.markdown("### 2D Structure")
            img = draw_molecule_2d(smiles)
            st.image(img)

        # Display 3D structure in the second column
        with col2:
            # SMILESから分子構造を生成してレンダリング
            blk = makeblock(smiles)
            st.markdown("### 3D Structure")
            render_mol(blk)
        
        st.markdown("### SDF(Structure-Data File)")
        st.code(blk)

    except Exception as e:
        st.warning(f"Unable to generate structure: {e}")


def pubchem_search():
    """PubChem APIを使った単分子分析アプリ。"""
    compound_name = st.text_input("名前を入力", value="acetone", key="name_input")
    if st.button("検索"):
        try:
            results = fetch_pubchem_data(compound_name)
            st.markdown("### CID number")
            st.code(results["cid"])
            st.markdown("### canonical_smiles")
            st.code(results["canonical_smiles"])
            st.markdown("### inchi")
            st.code(results["inchi"])
            st.markdown("### inchikey")
            st.code(results["inchikey"])  # InChIKeyを表示

            try:
                st.markdown("### 2D Structure")
                img = draw_molecule_2d(results["canonical_smiles"])
                st.image(img)
                # SMILESから分子構造を生成してレンダリング
                st.markdown("### 3D Structure")
                blk = makeblock(results["canonical_smiles"])
                render_mol(blk)
                st.markdown("### SDF(Structure-Data File)")
                st.code(blk)

            except Exception as e:
                # ユーザー向けエラーメッセージ
                st.error("3次元構造の描画中にエラーが発生しました。構造が正しいか確認してください。")
                st.error(f"エラー内容: {e}")

            with st.expander(f"📜 {compound_name}のpubchempyで取得できる内容一覧"):
                for index, data in results["data_frame"].items():
                    st.markdown(f"#### {index}")
                    st.code("\n".join(map(str, data.tolist())))
        except Exception as e:
            st.error(f"エラーが発生しました: {e}")

def smiles_to_data_display():
    # Streamlit app
    st.title("🔬 Display Structure and Molecular Properties from SMILES")

    # Example placeholder
    example_smiles = "CCO\nCC(=O)O\nC1=CC=CC=C1"

    # User input
    smiles_input = st.text_area("Paste SMILES here (one per line)", height=200, value=example_smiles)
    if st.button("Analyze"):
        # Process input
        smiles_list = [line.strip() for line in smiles_input.splitlines() if line.strip()]
        if smiles_list:
            st.info(f"Analyzing {len(smiles_list)} SMILES.")
            # Generate data
            data = smiles_to_data(smiles_list)

            # Display header row
            st.write("### Results")
            header_cols = st.columns([1, 2, 3, 2, 2])
            header_cols[0].write("**#**")
            header_cols[1].write("**SMILES**")
            header_cols[2].write("**Structure**")
            header_cols[3].write("**Molecular Weight**")
            header_cols[4].write("**molLogP**")

            # Display data in table format
            for index, entry in enumerate(data, start=1):
                col1, col2, col3, col4, col5 = st.columns([1, 2, 3, 2, 2])
                col1.write(f"#{index}")
                col2.write(entry["SMILES"])
                if isinstance(entry["Structure"], str):
                    col3.write(entry["Structure"])  # Error message if invalid
                else:
                    col3.image(entry["Structure"])  # Show structure image
                col4.write(f"{entry['MolecularWeight']:.2f}" if isinstance(entry["MolecularWeight"], float) else entry["MolecularWeight"])
                col5.write(f"{entry['molLogP']:.2f}" if isinstance(entry["molLogP"], float) else entry["molLogP"])

            # CSV download feature
            st.write("### Download Data")
            df = pd.DataFrame(
                [{"SMILES": d["SMILES"], "MolecularWeight": d["MolecularWeight"], "molLogP": d["molLogP"]} for d in data]
            )
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download CSV", data=csv, file_name="smiles_analysis.csv", mime="text/csv")
        else:
            st.warning("Please enter valid SMILES.")


if __name__ == "__main__":
    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()