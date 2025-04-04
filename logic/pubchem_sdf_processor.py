import os
from rdkit import Chem
import pandas as pd

def sdf_to_dataframe(sdf_file, selected_column=None):
    """
    SDFファイルを読み込み、指定した重要なカラムのみを含むPandasのDataFrameに変換する関数。
    無効な分子（SMILESが生成できない分子）はスキップします。

    Parameters:
        sdf_file (str): 入力するSDFファイルのパス。
        selected_column (list or str): 出力DataFrameに含めるカラムのリスト。
                                       "All"を指定すると全てのカラムを出力します。
                                       デフォルトは以下の4つ:
                                       ["PUBCHEM_CID_ASSOCIATIONS", "PUBCHEM_SUBSTANCE_SYNONYM",
                                        "PUBCHEM_EXT_SUBSTANCE_URL", "SMILES"].

    Returns:
        pd.DataFrame: 指定されたカラムのみを含む、または全てのカラムを含むSDFから抽出されたデータのDataFrame。
    """
    # デフォルトのカラムを設定
    if selected_column is None:
        selected_column = [
            "PUBCHEM_CID_ASSOCIATIONS",
            "PUBCHEM_SUBSTANCE_SYNONYM",
            "PUBCHEM_EXT_DATASOURCE_REGID",
            "PUBCHEM_EXT_SUBSTANCE_URL",
            "SMILES",
        ]

    # SDFファイルを読み込む
    supplier = Chem.SDMolSupplier(sdf_file)

    # データを格納するリストを作成
    data = []

    # 各分子の情報を抽出
    for mol in supplier:
        if mol is None:
            continue  # 無効な分子をスキップ
        
        # SMILESを生成
        smiles = Chem.MolToSmiles(mol)
        if not smiles:  # SMILESが空の場合はスキップ
            continue
        
        # 分子情報を取得
        mol_dict = mol.GetPropsAsDict()  # SDF中のプロパティを辞書形式で取得
        mol_dict["SMILES"] = smiles  # SMILESを追加
        data.append(mol_dict)

    # データをPandasのDataFrameに変換
    df = pd.DataFrame(data)

    # カラムの選択処理
    if selected_column != "All":  # "All"ではない場合、指定カラムのみ選択
        missing_columns = [col for col in selected_column if col not in df.columns]
        if missing_columns:
            print(f"以下のカラムは存在しません: {missing_columns}")
        selected_columns_in_df = [col for col in selected_column if col in df.columns]
        df = df[selected_columns_in_df]

    return df

def save_dataframe_in_chunks(df, base_filename, chunk_size_mb=50):
    """
    DataFrameを指定したサイズのチャンクに分割してCSVファイルとして保存する関数。

    Parameters:
        df (pd.DataFrame): 保存するDataFrame。
        base_filename (str): 基本となるCSVファイル名（例: 'output.csv'）。
        chunk_size_mb (int): 各CSVファイルの最大サイズ（MB単位）。
    """
    # 1行あたりの平均サイズを計算（概算）
    row_size = df.memory_usage(index=False, deep=True).sum() / len(df)
    max_rows = int((chunk_size_mb * 1024 * 1024) / row_size)

    # 分割して保存
    for i, chunk in enumerate(range(0, len(df), max_rows)):
        chunk_df = df.iloc[chunk:chunk + max_rows]
        chunk_filename = f"{os.path.splitext(base_filename)[0]}_part{i+1}.csv"
        chunk_df.to_csv(chunk_filename, index=False)
        print(f"CSVファイルが作成されました: {chunk_filename}")

