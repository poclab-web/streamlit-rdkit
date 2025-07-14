
# Chemoinformatics Tool 🧪

A comprehensive chemoinformatics learning and practice tool built with Streamlit and RDKit. Learn and explore a wide range of chemoinformatics methods from chemical structure analysis to machine learning.

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/streamlit-latest-red.svg)
![RDKit](https://img.shields.io/badge/rdkit-latest-green.svg)

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Key Libraries](#key-libraries)
- [Contributing](#contributing)
- [License](#license)

## 🎯 Overview

This tool is a comprehensive web application for learning chemoinformatics concepts from basic principles to advanced machine learning methods. It provides researchers, students, and chemical engineers with hands-on experience in molecular data analysis and prediction techniques.

## ✨ Features

### 1. Chemical Information (ChemInfo)
- Basic molecular structure editing and analysis
- SMILES notation processing
- Basic molecular property calculations

### 2. Structure Representation
- Convert molecular structures to various representation formats
- 2D/3D molecular structure visualization
- Interconversion between structure formats

### 3. Data Organization
- Chemical data cleaning and preprocessing
- Dataset integration and organization
- Data quality assessment

### 4. Chemical Numerical Values
- Physicochemical property calculations
- Statistical analysis of experimental values
- Numerical data visualization

### 5. Descriptors
- Molecular descriptor calculations
- Feature selection and evaluation
- Descriptor correlation analysis

### 6. Computational Chemistry
- Molecular orbital calculations
- Structure optimization
- Theoretical property calculations

### 7. Structure Search
- Substructure search
- Similarity search
- Structure pattern matching

### 8. Exploratory Analysis
- Data visualization
- Statistical analysis
- Pattern discovery

### 9. Evaluation Metrics
- Prediction accuracy evaluation
- Machine learning model performance assessment
- Regression and classification metrics

### 10. Machine Learning
- Molecular property prediction models
- Structure-activity relationship analysis
- Deep learning method applications

### 11. Interpretability
- Model interpretation
- Feature importance analysis
- SHAP value explanations

### 12. Difficulty Adjustment
- Content adjustment based on learning level
- Progressive learning support

### 13. Reverse Analysis
- Structure design from target properties
- Inverse problem solving methods

### 14. Chemical Space
- Chemical space visualization
- Diversity analysis
- Compound library evaluation

### 15. Interpolation/Extrapolation
- Data interpolation methods
- Extrapolation performance evaluation
- Applicability domain assessment

### 16. Literature Search
- Scientific literature search and analysis
- Research trend analysis

## 🚀 Installation

### Prerequisites
- Python 3.8 or higher
- pip or conda

### Installing Dependencies

```bash
# Clone the repository
git clone https://github.com/your-repo/streamlit-rdkit.git
cd streamlit-rdkit

# Install dependencies
pip install -r requirements.txt

# Or using conda
conda install -c conda-forge rdkit
pip install -r requirements.txt
```

### System Packages (Linux/macOS)
Some features require additional system packages:

```bash
# Based on packages.txt content
# Open Babel related dependencies
```

## 🖥️ Usage

### Starting the Application
```bash
streamlit run ChemoinformaticsTool.py
```

Access the application in your browser at `http://localhost:8501`.

### Basic Usage
1. Select the category you want to learn from the sidebar
2. Follow the instructions on each page
3. Experiment with sample data or your own data
4. Review and download results

## 📁 Project Structure

```
streamlit-rdkit/
├── ChemoinformaticsTool.py    # Main application
├── app_definitions.yaml       # App configuration and category definitions
├── constant.py               # Constants definition
├── requirements.txt          # Python dependencies
├── packages.txt             # System packages
├── README.md               # This file
├── common/                 # Common utilities
│   ├── display_app_overview.py
│   ├── multiapp.py
│   ├── sascore.py
│   └── stmolblock.py
├── data/                   # Data files
│   ├── curated-solubility-dataset.tab
│   ├── descriptors_name.csv
│   ├── FragranceDB.smi
│   ├── NMRshiftDB2_CHOonly_no_missing.xlsx
│   ├── PubChem_substance_phenol.sdf
│   ├── smiles.csv
│   ├── soac.csv
│   ├── TasteDB.smi
│   ├── TCI_fr.csv
│   ├── TCI_smiles.csv
│   └── Reagents/           # Reagent data
├── logic/                  # Business logic
│   ├── api.py
│   ├── chem_reactions.py
│   ├── chemical_search.py
│   ├── ChemRxiv_api.py
│   ├── FingerPrint.py
│   ├── IR_calculation.py
│   ├── load_data.py
│   ├── metrics.py
│   ├── mol_loader.py
│   ├── molecularconverter.py
│   ├── openbabel_utils.py
│   ├── pubchem_logic.py
│   ├── rdkit_draw_logic.py
│   ├── similarity.py
│   └── ...
├── pages/                  # Streamlit pages
│   ├── 01_ChemInfo.py
│   ├── 02_StructureRepresentation.py
│   ├── 03_DataOrganization.py
│   ├── 04_ChemicalNumericalValues.py
│   ├── 05_Descriptors.py
│   ├── 06_ComputationalChemistry.py
│   ├── 07_StructureSearch.py
│   ├── 08_ExploratoryAnalysis.py
│   ├── 09_EvaluationMetrics.py
│   ├── 10_MachineLearning.py
│   ├── 11_Interpretability.py
│   ├── 12_DifficultyAdjustment.py
│   ├── 13_ReverseAnalysis.py
│   ├── 14_ChemicalSpace.py
│   ├── 15_InterpolationExtrapolation.py
│   └── 16_LiteratureSearc.py
├── training/               # Machine learning model training
│   ├── gse.py
│   └── train_rf_models.py
├── utils/                  # Utilities
│   ├── app_runner.py
│   ├── display_code.py
│   ├── sidebar.py
│   ├── tab_handler.py
│   └── yaml_loader.py
└── vib/                   # Vibrational analysis
```

## 📚 Key Libraries

- **Streamlit**: Web application framework
- **RDKit**: Cheminformatics library
- **Pandas**: Data analysis
- **NumPy**: Numerical computing
- **Matplotlib/Seaborn**: Data visualization
- **Scikit-learn**: Machine learning
- **PyCaret**: AutoML
- **SHAP**: Model interpretation
- **Py3Dmol**: 3D molecular visualization
- **PubChemPy**: PubChem API
- **OpenBabel**: Chemical format conversion
- **SELFIES**: Molecular representation

## 🤝 Contributing

### How to Contribute
1. Fork this repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Create a Pull Request

### Development Guidelines
- Add new pages to the `pages/` directory with modular design
- Place common functions in `common/` or `utils/`
- Organize business logic in `logic/`
- Add new dependencies to `requirements.txt`
- Include appropriate comments in code
- Design for testability

### Bug Reports & Feature Requests
Please follow the Issue templates when reporting.

## 📄 License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.

## 🌟 Acknowledgments

- RDKit development team
- Streamlit development team
- Cheminformatics community

## 📞 Support

If you have questions or issues, please contact us through:
- GitHub Issues
- Pull Requests
- Discussions

---

⚠️ **Note**: This application is currently under development. We aim to complete version 1.0 by summer 2025.
