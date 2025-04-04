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
    # Streamlit アプリ
    st.title("🔬 SMILESから構造と分子特性を表示")

    # プレースホルダーに例を設定
    example_smiles = "CCO\nCC(=O)O\nC1=CC=CC=C1"

    # ユーザー入力
    smiles_input = st.text_area("SMILESを貼り付けてください（複数行）", height=200, value=example_smiles)
    if st.button("解析"):
        # 入力を処理
        smiles_list = [line.strip() for line in smiles_input.splitlines() if line.strip()]
        if smiles_list:
            st.info(f"{len(smiles_list)} 件のSMILESを解析しています。")
            # データを生成
            data = smiles_to_data(smiles_list)

            # ヘッダー行を表示
            st.write("### 結果一覧")
            header_cols = st.columns([1, 2, 3, 2, 2])
            header_cols[0].write("**#**")
            header_cols[1].write("**SMILES**")
            header_cols[2].write("**構造**")
            header_cols[3].write("**分子量**")
            header_cols[4].write("**molLogP**")

            # データをテーブル形式で表示
            for index, entry in enumerate(data, start=1):  # 1から始まるインデックス
                col1, col2, col3, col4, col5 = st.columns([1, 2, 3, 2, 2])
                col1.write(f"#{index}")  # インデックス番号を表示
                col2.write(entry["SMILES"])
                if isinstance(entry["構造"], str):
                    col3.write(entry["構造"])  # 無効な場合はエラーメッセージ
                else:
                    col3.image(entry["構造"])  # 構造画像を表示
                col4.write(f"{entry['分子量']:.2f}" if isinstance(entry["分子量"], float) else entry["分子量"])
                col5.write(f"{entry['molLogP']:.2f}" if isinstance(entry["molLogP"], float) else entry["molLogP"])

            # CSVデータのダウンロード機能
            st.write("### データをダウンロード")
            df = pd.DataFrame(
                [{"SMILES": d["SMILES"], "MolWt": d["分子量"], "molLogP": d["molLogP"]} for d in data]
            )
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 CSVをダウンロード", data=csv, file_name="smiles_analysis.csv", mime="text/csv")
        else:
            st.warning("有効なSMILESを入力してください。")


if __name__ == "__main__":
    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()