import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi
import joblib
from sklearn.metrics import r2_score, mean_absolute_error
import matplotlib.pyplot as plt

# アプリの定義

@st.cache_data
def load_shift_output():
    return pd.read_excel(".models/shift_output_file.xlsx")

def nmr_13c_hose_display():
    st.title("📊 NMRのデータとHOSEコード")
    st.write("ここではNMRのデータとHOSEコードを表示します。")
    df = load_shift_output()

    smiles_list = df["smiles"].unique()
    num_compounds = len(smiles_list)
    st.info(f"化合物の種類: {num_compounds} 種類")

    idx = st.number_input("表示する化合物の番号を選んでください（1〜{0}）".format(num_compounds), min_value=1, max_value=num_compounds, value=1)
    smiles = smiles_list[idx-1]

    st.markdown(f"## 分子 {idx}: `{smiles}`")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # 分子式
        formula = rdMolDescriptors.CalcMolFormula(mol)
        st.write(f"**分子式:** {formula}")

        # InChIKey（示性式の代用例）
        inchikey = inchi.MolToInchiKey(mol)
        st.write(f"**InChIKey:** {inchikey}")

        mol = Chem.MolFromSmiles(smiles)
        st.image(Draw.MolToImage(mol, size=(300, 300)), caption=smiles)


    display_cols = [
        "Shift_ID", "Shift_Value", "Atom_Index", "Atom_Type",
        "HOSE_Code_R1", "HOSE_Code_R2", "HOSE_Code_R3"
    ]
    group = df[df["smiles"] == smiles]
    st.dataframe(group[display_cols].reset_index(drop=True))


@st.cache_data
def load_cho_compounds():
    return pd.read_csv(".models/filtered_cho_compounds.csv", encoding='utf-8')

def structure_isomer_output_display():
    st.title("📊 構造異性体の出力")
    st.write("ここでは、構造異性体の出力を行います。")
    df = load_cho_compounds()

    # SMILES列から分子式を計算して新しい列を追加
    if "smiles" in df.columns:
        df["formula"] = df["smiles"].apply(lambda x: rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else "")
    else:
        st.error("SMILES列が見つかりません。")
        return

    # ユーザーに分子式を入力させる
    formula_input = st.text_input("分子式を入力してください（例: C2H6O）")

    import re

    def parse_formula(formula):
        # 元素記号と数を抽出
        matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        elements = {}
        for elem, num in matches:
            elements[elem] = int(num) if num else 1
        return elements

    if formula_input:
        elements = parse_formula(formula_input)
        allowed_elements = {"C", "H", "O"}
        # C, H, O以外が含まれている場合
        if not set(elements.keys()).issubset(allowed_elements):
            st.warning("C、H、O以外の元素が含まれています。このデータには対応していません。")
        # CとOの合計が9を超える場合
        elif elements.get("C", 0) + elements.get("O", 0) > 9:
            st.warning("CとOの合計が9を超えています。このデータには対応していません。")
        else:
            filtered_df = df[df["formula"] == formula_input]
            if filtered_df.empty:
                st.info(f"分子式 `{formula_input}` の構造異性体はデータに存在しません。")
            else:
                st.write(f"分子式 `{formula_input}` の構造異性体: {len(filtered_df)} 件")
                for idx, row in filtered_df.iterrows():
                    st.markdown(f"### SMILES: `{row['smiles']}`")
                    mol = Chem.MolFromSmiles(row["smiles"])
                    if mol:
                        st.image(Draw.MolToImage(mol, size=(200, 200)), caption=row["smiles"])
                    st.write(row)
    else:
        st.info("分子式を入力してください。")


# TODO: 逆解析の表示機能を実装する
def reverse_analysis_display():
    st.title("📊 逆解析")
    st.write("ここでは、予測モデルを組み合わせて構造を出力します。")

    st.warning("この機能はまだ実装されていません。")

# 以下は実装例のコメントアウト
#     # データの読み込み
#     df = load_shift_output()
#     if "smiles" not in df.columns or "Shift_Value" not in df.columns:
#         st.error("必要なカラム（smiles, Shift_Value）がデータにありません。")
#         return

#     # 特徴量生成（例：分子記述子を使う場合）
#     from rdkit.Chem import Descriptors

#     def featurize(smiles):
#         mol = Chem.MolFromSmiles(smiles)
#         if mol:
#             return [
#                 Descriptors.MolWt(mol),
#                 Descriptors.NumHDonors(mol),
#                 Descriptors.NumHAcceptors(mol),
#                 Descriptors.MolLogP(mol)
#             ]
#         else:
#             return [0, 0, 0, 0]

#     X = df["smiles"].apply(featurize).tolist()
#     y = df["Shift_Value"].values

#     # モデルの読み込み
#     model = joblib.load(".models/random_forest_model.pkl")

#     # 予測
#     y_pred = model.predict(X)

#     # 精度指標
#     r2 = r2_score(y, y_pred)
#     mae = mean_absolute_error(y, y_pred)
#     st.write(f"**R²スコア:** {r2:.3f}")
#     st.write(f"**MAE:** {mae:.3f}")

#     # yyplot
#     fig, ax = plt.subplots()
#     ax.scatter(y, y_pred, alpha=0.7)
#     ax.plot([y.min(), y.max()], [y.min(), y.max()], 'r--')
#     ax.set_xlabel("実測値")
#     ax.set_ylabel("予測値")
#     ax.set_title("yyplot（予測値 vs. 実測値）")
#     st.pyplot(fig)


if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "ReverseAnalysis"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()