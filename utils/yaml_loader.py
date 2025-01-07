import yaml

def load_yaml(file_path: str) -> dict:
    """
    YAMLファイルを読み込み、辞書形式で返す関数。

    Parameters:
        file_path (str): YAMLファイルのパス。

    Returns:
        dict: YAMLファイルの内容を辞書形式で返す。
    """
    with open(file_path, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)
