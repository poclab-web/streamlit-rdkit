import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from streamlit_ketcher import st_ketcher

import pandas as pd
import io
from logic.rdkit_salt_stereo_analyzer import MoleculeAnalyzer, analyze_chirality
from logic.chem_reactions import ReactionSmilesParser

from rdkit.Chem import Draw
from rdkit import Chem

# アプリの定義

def molecule_analyzer_display():
    # Streamlit App
    st.title("Molecule Analyzer")

    # アプリ調整中の警告文
    st.warning("このアプリ(Molecule Analyzer)は現在調整中です。一部の機能(立体判定)が正しく動作しない可能性があります。")

    # 複数のSMILES入力フィールド
    smiles_input = st.text_area("Enter one or more SMILES strings (one per line) for analysis")

    # SMILES解析ボタン
    analyze_smiles = st.button("Analyze SMILES")

    # ファイルアップロード
    uploaded_file = st.file_uploader("Upload an SDF or CSV file", type=["sdf", "csv"])

    # デフォルトデモファイルのパス
    demo_file_path = "data/PubChem_substance_phenol.sdf"

    # デモ実行用のボタン
    use_demo_data = st.button("Run Analysis with Demo Data")

    # SMILES解析
    if analyze_smiles and smiles_input.strip():
        with st.spinner("Analyzing SMILES strings..."):
            try:
                smiles_list = [smiles.strip() for smiles in smiles_input.split("\n") if smiles.strip()]
                data = []
                for smiles in smiles_list:
                    analyzer = MoleculeAnalyzer(source_type='smiles', source=smiles)
                    df, _ = analyzer.analyze()
                    data.append(df)
                # 複数の結果を結合
                result_df = pd.concat(data, ignore_index=True)
                summary = {
                    "Total Molecules": len(result_df),
                    "Molecules with Salts": result_df['HasSalt'].sum(),
                    "Molecules with Undefined Stereo": result_df['HasUndefinedStereo'].sum(),
                    "Unique Molecules": len(result_df['SMILES'].unique()),
                    "Duplicate Molecules": result_df['IsDuplicate'].sum()
                }
            except Exception as e:
                st.error(f"Error analyzing SMILES: {e}")
                return

    # ファイルアップロード解析
    elif uploaded_file is not None:
        if uploaded_file.name.endswith(".sdf"):
            with st.spinner("Analyzing SDF file..."):
                file_obj = io.BytesIO(uploaded_file.read())
                analyzer = MoleculeAnalyzer(source_type='sdf', source=file_obj)
                result_df, summary = analyzer.analyze()
        elif uploaded_file.name.endswith(".csv"):
            with st.spinner("Analyzing CSV file..."):
                csv_df = pd.read_csv(uploaded_file)
                analyzer = MoleculeAnalyzer(source_type='csv', source=csv_df)
                result_df, summary = analyzer.analyze()

    # デモデータ解析
    elif use_demo_data:
        st.info(f"Using demo data from '{demo_file_path}'.")
        with open(demo_file_path, "rb") as demo_file:
            with st.spinner("Analyzing demo SDF file..."):
                file_obj = io.BytesIO(demo_file.read())
                analyzer = MoleculeAnalyzer(source_type='sdf', source=file_obj)
                result_df, summary = analyzer.analyze()

    else:
        st.warning("Please enter SMILES strings, upload a file, or click the button to use the demo data.")
        return

    # 解析結果の表示
    st.subheader("Analysis Summary")
    for key, value in summary.items():
        st.write(f"{key}: {value}")

    # 脱塩や脱立体化後の重複を表示
    st.subheader("Duplicates After Desalting and Stereochemistry Removal")
    try:
        # 脱塩および脱立体化したSMILESを生成
        result_df['Desalted_SMILES'] = result_df['SMILES'].apply(lambda x: MoleculeAnalyzer.desalt_smiles(x))
        result_df['Canonical_SMILES'] = result_df['Desalted_SMILES'].apply(lambda x: MoleculeAnalyzer.remove_stereo(x))

        # 重複の確認
        unique_canonical_smiles = result_df['Canonical_SMILES'].nunique()
        duplicate_count = len(result_df) - unique_canonical_smiles

        st.write(f"Unique Molecules After Processing: {unique_canonical_smiles}")
        st.write(f"Duplicate Molecules After Processing: {duplicate_count}")

        # 脱塩および脱立体化後のデータフレームを表示
        st.subheader("Processed Molecule Data")
        st.dataframe(result_df[['SMILES', 'Desalted_SMILES', 'Canonical_SMILES', 'HasSalt', 'HasUndefinedStereo']])
    except Exception as e:
        st.error(f"Error processing duplicates after desalting and stereochemistry removal: {e}")

    # 立体化学情報のデータフレーム表示
    st.subheader("Stereochemistry Details")
    try:
        # 立体化学情報を抽出
        stereo_data = []
        for _, row in result_df.iterrows():
            undefined_stereo = row.get("UndefinedStereo", {})
            chiral_centers = undefined_stereo.get("undefined_chiral_centers", [])
            ez_bonds = undefined_stereo.get("undefined_ez_bonds", [])

            # 不斉点情報をデータフレーム用に整形
            for center in chiral_centers:
                stereo_data.append({
                    "SMILES": row["SMILES"],
                    "Type": "Chiral Center",
                    "Index": center["index"],
                    "Details": center["type"]
                })

            # 二重結合情報をデータフレーム用に整形
            for bond in ez_bonds:
                stereo_data.append({
                    "SMILES": row["SMILES"],
                    "Type": "Double Bond",
                    "Index": f"{bond['begin_atom_idx']}-{bond['end_atom_idx']}",
                    "Details": "Undefined"
                })

        # データフレームに変換
        stereo_df = pd.DataFrame(stereo_data)

        # データフレームを表示
        if not stereo_df.empty:
            st.dataframe(stereo_df)
        else:
            st.info("No undefined stereochemistry found.")
    except Exception as e:
        st.error(f"Error displaying stereochemistry details: {e}")

    # DataFrame の表示
    st.subheader("Molecule Data")
    st.dataframe(result_df[['SMILES', 'HasSalt', 'HasUndefinedStereo', 'IsDuplicate']])

    # DataFrame をダウンロード可能に
    csv = result_df.to_csv(index=False)
    st.download_button(
        label="Download DataFrame as CSV",
        data=csv,
        file_name="analyzed_data.csv",
        mime="text/csv"
    )

def checkStereocenters():
    # Streamlitアプリの開始
    st.title("不斉点判定アプリ")
    st.write("構造を描画して不斉点を判定し、CIP順位を表示します。")

    # 構造をKetcherで入力
    smiles = st_ketcher()

    # ユーザーが構造を入力したか確認
    if not smiles:
        st.warning("構造が描画されていません。構造を描画してApplyをクリックしてください。")
    else:
        # 入力されたSMILESの表示
        st.write("入力されたSMILES:")
        st.code(smiles)

        # 不斉点の解析
        mol, chiral_centers, unassigned, chiral_details = analyze_chirality(smiles)

        if mol:
            # インデックス番号付きの分子構造を描画
            st.write("### 分子構造とインデックス")
            mol_with_indices = Chem.Mol(mol)
            for atom in mol_with_indices.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
            
            img = Draw.MolToImage(mol_with_indices, size=(400, 400))
            st.image(img, caption="インデックス番号付きの分子構造")

            # 不斉点が存在するかを確認
            if chiral_centers:
                st.success(f"不斉点が {len(chiral_centers)} 箇所見つかりました。")
                
                # 不斉点の詳細を表示
                for detail in chiral_details:
                    st.write(f"### 不斉点（原子インデックス {detail['center_idx']}）")
                    st.write(f"- 立体化学: {'未割り当て' if detail['label'] == '?' else detail['label']}")
                    st.write(f"- 置換基:")
                    for sub in detail["substituents"]:
                        st.write(f"  - 原子インデックス {sub[0]}: 元素 {sub[1]}, 質量 {sub[2]:.2f}")
                
                # 立体表記のない不斉点があるか
                if unassigned:
                    st.warning("立体表記がない不斉点があります。")
                    for idx, _ in unassigned:
                        st.write(f"原子インデックス {idx} の立体表記がありません。")
            else:
                st.info("この分子には不斉点がありません。")


        else:
            st.error("分子構造の解析に失敗しました。")

def reaction_analyzer_display():
    st.title("REACTION SMILES Parser")
    st.write("REACTION SMILES を解析して Reactants, Reagents, Products を表示します。")

    # デフォルトのREACTION SMILES
    default_smiles = """CCO.CCOC>>O=CC.CCO
CCO>>CC=O.CC
CCC.OO>O=O>CC(C)O"""

    # ユーザー入力
    smiles_input = st.text_area("REACTION SMILESを入力してください (複数行入力可能):", default_smiles)

    if st.button("解析"):  # ボタンを押すと解析を実行
        if smiles_input.strip():
            # 入力を改行で分割してリスト化
            reaction_smiles_list = smiles_input.strip().split("\n")

            # ReactionSmilesParser を使用して解析
            parser = ReactionSmilesParser(reaction_smiles_list)
            df = parser.to_dataframe()

            # データフレームを表示
            st.write("### 解析結果")
            st.dataframe(df)

            # Reactants, Reagents, Productsを描画して表示
            st.write("### 分子の可視化")
            for index, row in df.iterrows():
                cols = st.columns(3)
                with cols[0]:
                    st.write(f"**Reactants ({index + 1})**")
                    reactant_img = parser.draw_molecules(row['Reactants'])
                    st.image(reactant_img, caption="Reactants", use_container_width=False)

                with cols[1]:
                    if row['Reagents']:
                        st.write(f"**Reagents ({index + 1})**")
                        reagent_img = parser.draw_molecules(row['Reagents'])
                        st.image(reagent_img, caption="Reagents", use_container_width=False)

                with cols[2]:
                    st.write(f"**Products ({index + 1})**")
                    product_img = parser.draw_molecules(row['Products'])
                    st.image(product_img, caption="Products", use_container_width=False)
        else:
            st.error("REACTION SMILESを入力してください！")

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "DataOrganization"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()