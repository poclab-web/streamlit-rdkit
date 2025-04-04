# ロジック部分
import requests

def fetch_chemrxiv_data(keyword):
    """
    ChemRxiv APIから指定されたキーワードに基づくデータを取得します。

    Parameters:
        keyword (str): 検索キーワード。

    Returns:
        dict: APIのレスポンスデータ。
    """
    url = "https://chemrxiv.org/engage/chemrxiv/public-api/v1/items"
    params = {"term": f'"{keyword}"'}

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}

def parse_chemrxiv_data(data):
    """
    ChemRxiv APIから取得したデータを解析してフォーマットします。

    Parameters:
        data (dict): ChemRxiv APIのレスポンスデータ。

    Returns:
        dict: フォーマット済みのデータ。
    """
    if "error" in data:
        return data

    parsed_data = {
        "totalCount": data.get("totalCount", 0),
        "items": []
    }

    for item_hit in data.get("itemHits", []):
        item = item_hit.get("item", {})
        parsed_data["items"].append({
            "title": item.get("title", "タイトルなし"),
            "doi": item.get("doi", "DOIなし"),
            "authors": ", ".join(
                f"{author.get('firstName', '')} {author.get('lastName', '')}"
                for author in item.get("authors", [])
            ),
            "publishedDate": item.get("publishedDate", "公開日不明"),
            "abstract": item.get("abstract", "概要なし"),
            "keywords": ", ".join(item.get("keywords", [])),
            "license_name": item.get("license", {}).get("name", "ライセンス情報なし"),
            "license_url": item.get("license", {}).get("url", "URLなし")
        })

    return parsed_data
