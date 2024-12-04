import streamlit as st
from pycaret.regression import *
import constant

def search_data_types(setup_result, list_columns: list, target: str):
    """
    Search input of data types
    Args:
        setup_result: tuple
            transformation pipeline
        list_columns: list
            list of columns in dataset
        target: str
            target column
    Returns:
        list_num: list
            list of strings with column names that are numeric.
        list_cat: list
            list of strings with column names that are categorical.
    """

    list_num_org = list_columns.copy()

    for i in range(len(setup_result)):
        try:
            _ = list(setup_result[i].columns)
        except:
            pass
        else:
            if list(setup_result[i].columns) != list_columns:
                list_num_org = list(setup_result[i].columns)
                break
            else:
                pass

    list_num = list(set(list_num_org) & set(list_columns))

    list_cat = list_columns.copy()
    for i in list_num:
        list_cat.remove(i)

    list_cat.remove(target)

    return list_num, list_cat


@st.cache_data(allow_output_mutation=True)
def compare_regression_models():
    """
    Run compare_models
    Args: none
    Returns: pd.DataFrame
        result
    """
    best = compare_models()
    best_model_results = pull()
    return pd.DataFrame(best_model_results).reset_index(drop=True)


@st.cache_data(allow_output_mutation=True)
def create_regression_model(estimator):
    return create_model(estimator)


@st.cache_data(allow_output_mutation=True)
def create_regression_blend_model(list_estimator):
    _list_models = []
    for _model in list_estimator:
        _list_models.append(create_model(_model))
    return blend_models(estimator_list=_list_models)


def create_model_sequence(data):
    """
    create model
    Args:
        data: pd.DataFrame
            dataset
    Returns: model
    """

    model = None
    list_columns = list(data.columns)

    # Select some options for the setup function.
    with st.form("select_options"):
        # required
        target_column = st.selectbox("target variable", list_columns)

        # optional
        with st.expander("Options"):
            done_normalize = st.checkbox("done normalize")
            done_pca = st.checkbox("done dimensionality reduction")

            num_train_size = st.number_input(
                "train size",
                min_value=0.0,
                max_value=1.0,
                value=0.7,
                step=0.05,
                help="Proportion of the dataset to be used for training and validation.",
            )

            list_cat = st.multiselect("categorical features", list_columns)
            list_num = st.multiselect("numeric features", list_columns)
            list_igr = st.multiselect("ignore features", list_columns)

        done_pycaret = st.form_submit_button("done")

        if done_pycaret:
            st.session_state.done_pycaret = True

    if st.session_state.done_pycaret:
        # setup
        pipeline = setup(
            data,
            target=target_column,
            train_size=num_train_size,
            categorical_features=list_cat,
            numeric_features=list_num,
            ignore_features=list_igr,
            normalize=done_normalize,
            pca=done_pca,
            html=False,  # Streamlit
            silent=True,  # Streamlit
        )

        # confirmation input of data types
        list_num, list_cat = search_data_types(pipeline, list_columns, target_column)
        st.write("Numeric Features:")
        st.info("%s" % str(list_num)[1:-1])
        st.write("Categorical Features:")
        st.info("%s" % str(list_cat)[1:-1])

        # compare models
        st.markdown("---")
        st.subheader("Compare Models")
        df_result_compare = compare_regression_models()
        st.write("Result of compare models")
        st.write(df_result_compare)

        list_choise_name = []
        list_choise_id = []

        with st.form("select_using_models"):
            list_choise_name = st.multiselect(
                "Choose any models",
                list(df_result_compare["Model"]),
                help="If multiple models are selected, a blended model will be generated.",
            )

            done_choise = st.form_submit_button("done")

        if done_choise:
            # Convert model names and IDs
            for _model_name in list_choise_name:
                list_choise_id.append(constant.REGRESSION_MODELS[_model_name])

            if len(list_choise_id) == 0:  # NONE
                st.warning("Please select some models.")
            elif len(list_choise_id) == 1:  # single
                st.session_state["model"] = create_regression_model(list_choise_id[0])
            else:  # blendmodel
                st.session_state["model"] = create_regression_blend_model(list_choise_id)

        else:  # no select models
            pass  # DO NOTHING

    else:  # no select setup options
        pass  # DO NOTHING


@st.cache_data
def convert_df(df):
   return df.to_csv().encode('utf-8')


def pycaret_regression():
    # Setting Streamlit
    # st.set_page_config(page_title="PyCaret Demo", layout="wide")

    st.title('Pycaret 😀')

    # load dateset
    uploaded_file = st.file_uploader("csvファイルをアップロードしてください")

    # test dataset
    test_df = pd.read_csv("data/soac.csv")
    test = convert_df(test_df)

    st.download_button(
        "example csvのDownload",
        test,
        "example.csv",
        "text/csv",
        key='download-csv'
    )

    if uploaded_file is not None:
        uploaded_data = pd.read_csv(uploaded_file)
        st.write(uploaded_data)
        # create model
        create_model_sequence(uploaded_data)

