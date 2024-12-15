import pubchempy as pcp
import pandas as pd

def fetch_pubchem_data(compound_name="acetone"):
    """
    PubChemからデータ("canonical_smiles", "cid", DataFrame)を取得するロジック部分

    Requirements:
    この関数を使うためには、pubchempyが必要です。
        
    ```
    pip install pubchempy
    ```

    Parameters:
        compound_name (str): 検索する化合物の名前 (CAS番号、慣用名、IUPAC名など)

    Returns:
        dict: 検索結果を含む辞書形式のデータ
        例: {
            "canonical_smiles": "CCO",
            "cid": 12345,
            "data_frame": DataFrame
        }

    Raises:
        ValueError: データが存在しない場合
        Exception: その他のエラー
    
    """
    try:
        # PubChemからデータを取得
        df = pcp.get_compounds(compound_name, 'name', as_dataframe=True)

        # データが存在するかを確認
        if df.empty:
            raise ValueError("該当するデータがありません。")

        # 必要なデータを抽出
        canonical_smiles = df.at[df.index.values[0], "canonical_smiles"]
        cid = df.index.values[0]

        return {
            "canonical_smiles": canonical_smiles,
            "cid": cid,
            "data_frame": df
        }
    except ValueError as ve:
        raise ve
    except Exception as e:
        raise Exception(f"PubChemデータ取得中にエラーが発生しました: {e}")
