import streamlit as st
from utils.tab_handler import handle_tabs_for_category
from utils.sidebar import display_sidebar

from logic.rdkit_salt_stereo_analyzer import MoleculeAnalyzer


# アプリの定義

def molecule_analyzer_display():
    # Streamlit App
    st.title("Molecule Analyzer")

    uploaded_file = st.file_uploader("Upload an SDF or CSV file", type=["sdf", "csv"])

    if uploaded_file is not None:
        if uploaded_file.name.endswith(".sdf"):
            with st.spinner("Analyzing SDF file..."):
                analyzer = MoleculeAnalyzer(source_type='sdf', source=uploaded_file)
                df, summary = analyzer.analyze()
        elif uploaded_file.name.endswith(".csv"):
            with st.spinner("Analyzing CSV file..."):
                csv_df = pd.read_csv(uploaded_file)
                analyzer = MoleculeAnalyzer(source_type='csv', source=csv_df)
                df, summary = analyzer.analyze()

        # Display summary
        st.subheader("Analysis Summary")
        for key, value in summary.items():
            st.write(f"{key}: {value}")

        # Display DataFrame
        st.subheader("Molecule Data")
        st.dataframe(df[['SMILES', 'HasSalt', 'HasUndefinedStereo', 'IsDuplicate']])

        # Allow download of DataFrame
        csv = df.to_csv(index=False)
        st.download_button(
            label="Download DataFrame as CSV",
            data=csv,
            file_name="analyzed_data.csv",
            mime="text/csv"
        )

if __name__ == "__main__":
    # 現在のカテゴリー（手動設定）
    current_category = "DataOrganization"  # 正しいカテゴリーキーを指定
    st.write(f"現在のカテゴリー: {current_category}")  # デバッグ用

    # ページ共通のタブ処理
    handle_tabs_for_category(current_category)

    # サイドバーを表示
    display_sidebar()