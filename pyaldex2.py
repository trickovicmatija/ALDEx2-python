import pandas as pd

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter



def run_aldex2(counts:pd.DataFrame,metadata:pd.DataFrame,test:str,r_script_path:str,mc_samples='auto')->pd.DataFrame:
    """
    Run ALDEx2 on a raw counts dataframe.
    
    Returns:
        pd.DataFrame: A dataframe with the results of ALDEx2.
    
    Arguments:
        counts -> pd.DataFrame: A dataframe with the raw counts.
        metadata -> pd.DataFrame: A dataframe with the metadata.
        test -> str: The test to run. Can be one of the following:
            't': Welch's t-test together with non-parametric Wilcoxon test,
            'kw': Kruskal-Wallis test together with glm test,
            'glm': Generalized linear model test using model.matrix,
            'corr': Correlation test using cor.test.
        r_script_path -> str: The path to the R script to run ALDEx2.
        mc_samples -> int: The number of Monte Carlo samples to use. If set to 'auto',
            the number of samples is determined automatically.
 
        """
    r = ro.r
    r['source'](str(r_script_path))
    run_aldex_function_r = ro.globalenv['run_aldex']
    print("Successfully initiated R function!")
    if counts.shape[1] != metadata.shape[0]:
        counts = counts.T
    assert counts.shape[1] == metadata.shape[0], "Counts dataframe and metadata need to have same number of samples"
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_r = ro.conversion.py2rpy(counts)
    metadata = metadata.loc[counts.columns].iloc[:,0]
    categories = metadata.values.tolist()
    cond = ro.StrVector(categories)
    if mc_samples=='auto' :
        montecarlo_samples = int(1000 / metadata.value_counts().values.min()) # This is the lowest number of samples recommended by aldex2.
    else:
        montecarlo_samples = int(mc_samples)
    print(f"Running with {montecarlo_samples} montecarlo samples!")
    print("Starting ALDEx2!")
    df_result_r = run_aldex_function_r(df_r, cond, montecarlo_samples, test)
    print("Finished ALDEx2!")
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_result = ro.conversion.rpy2py(df_result_r)
    return df_result