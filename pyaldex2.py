import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
        test -> str: The test to run. Can be one of the following:s
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
    input_counts = counts.copy()
    if input_counts.shape[1] != metadata.shape[0]:
        input_counts = input_counts.T
    assert input_counts.shape[1] == metadata.shape[0], "Counts dataframe and metadata need to have same number of samples"
    metadata = metadata.loc[input_counts.columns].iloc[:,0]
    categories = metadata.values.tolist()
    cond = ro.StrVector(categories)
    if mc_samples=='auto' :
        montecarlo_samples = int(1000 / metadata.value_counts().values.min()) # This is the lowest number of samples recommended by aldex2.
    else:
        montecarlo_samples = int(mc_samples)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_r = ro.conversion.py2rpy(input_counts)
    print(f"Running with {montecarlo_samples} montecarlo samples!")
    print("Starting ALDEx2!")
    df_result_r = run_aldex_function_r(df_r, cond, montecarlo_samples, test)
    print("Finished ALDEx2!")
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_result = ro.conversion.rpy2py(df_result_r)
    return df_result


def MA_plot(results_input,effect_threshold=2,horizontal_line=True,figsize=(12,8),title='MA plot'):
    """
    Plot the MA plot of ALDEx2 results.

    Returns:
        None.

    Arguments:
        results_input -> pd.DataFrame: A dataframe with the results of ALDEx2.
        effect_threshold -> float: The threshold for the effect size value in order to mark the points red.
        horizontal_line -> bool: If True, a horizontal line is drawn at value 0.
        figsize -> tuple: The size of the figure.
        title -> str: The title of the figure.
    """
    results = results_input.copy()
    results['effect_hue']=results['effect'].apply(lambda x: f'abs(effect) > {effect_threshold}' if abs(x) >effect_threshold else f'abs(effect) < {effect_threshold}')
    results['size']=results['effect'].apply(lambda x: 'big' if abs(x) >effect_threshold else "small")
    color_dict = dict({f'abs(effect) > {effect_threshold}':'red',f'abs(effect) < {effect_threshold}': 'black'})
    size_dict=dict({'small':5,'big':20})
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(x='rab.all',y='diff.btw',hue='effect_hue',data=results,palette=color_dict,size='size',sizes=size_dict,legend='full',ax=ax)
    if horizontal_line:
        plt.axhline(y=0,color='black',lw=1, ls='--')
    legend = ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    if len(results['size'].value_counts()) == 1:
        ax.get_legend().remove()
    else:
        ax.legend(handles=handles[1:3], labels=labels[1:3], title="")
    plt.ylabel('Difference')
    plt.xlabel('Abundance')
    plt.title(title)
    return None