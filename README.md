




# 参考ファイル
[PyCaretとStreamlitでAutoMLのGUIツールをさくっと作ってみる](https://qiita.com/well_living/items/e9560ade1f0adedbaf6c)

[Streamlit上でPyCaretを動かす方法](https://qiita.com/nockn/items/77d6b5f5e8f58a0b6c44)

## Improvements and Guidelines

### Project Structure
- `main.py`: Entry point for the Streamlit application.
- `common/`: Utility functions and modules used across the app.
- `data/`: Contains data files such as CSV and SDF for analysis.
- `main_pages/`: Individual pages of the Streamlit multi-page app.

### Setup Instructions
1. Clone the repository.
2. Install dependencies: `pip install -r requirements.txt`.
3. Run the application: `streamlit run main.py`.

### Contribution Guidelines
- Ensure new pages are modular and added to the `main_pages/` directory.
- Keep utility functions in `common/`.
- Add new dependencies to `requirements.txt`.
