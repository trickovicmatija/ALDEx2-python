import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys


import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def run_aldex2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    test: str,
    r_script_path: str,
    mc_samples="auto",
) -> pd.DataFrame:
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
    r["source"](str(r_script_path))
    run_aldex_function_r = ro.globalenv["run_aldex"]
    print("Successfully initiated R function!")
    input_counts = counts.copy()
    if input_counts.shape[1] != metadata.shape[0]:
        input_counts = input_counts.T
    assert (
        input_counts.shape[1] == metadata.shape[0]
    ), "Counts dataframe and metadata need to have same number of samples"
    metadata = metadata.loc[input_counts.columns].iloc[:, 0]
    categories = metadata.values.tolist()
    cond = ro.StrVector(categories)
    if mc_samples == "auto":
        montecarlo_samples = int(
            1000 / metadata.value_counts().values.min()
        )  # This is the lowest number of samples recommended by aldex2.
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


def get_clr_instance(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    r_script_path: str,
    instance=1,
    mc_samples="auto",
) -> pd.DataFrame:
    """
    Run ALDEx2 to get one CLR instance.

    Returns:
        pd.DataFrame: A dataframe with the CLR-transformed data.

    Arguments:
        counts -> pd.DataFrame: A dataframe with the raw counts.
        metadata -> pd.DataFrame: A dataframe with the metadata.
        r_script_path -> str: The path to the R script to run ALDEx2.
        mc_samples -> int: The number of Monte Carlo samples to use. If set to 'auto',
            the number of samples is determined automatically.
    """
    r = ro.r
    r["source"](str(r_script_path))
    get_clr_function_r = ro.globalenv["get_clr"]
    input_counts = counts.copy()
    if input_counts.shape[1] != metadata.shape[0]:
        input_counts = input_counts.T
    assert (
        input_counts.shape[1] == metadata.shape[0]
    ), "Counts dataframe and metadata need to have same number of samples"
    metadata = metadata.loc[input_counts.columns].iloc[:, 0]
    categories = metadata.values.tolist()
    cond = ro.StrVector(categories)
    if mc_samples == "auto":
        montecarlo_samples = int(
            1000 / metadata.value_counts().values.min()
        )  # This is the lowest number of samples recommended by aldex2.
    else:
        montecarlo_samples = int(mc_samples)
    input_counts = input_counts.loc[
        input_counts.sum(1) > 0
    ]  # Remove features with zero counts in all samples.
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_r = ro.conversion.py2rpy(input_counts)
    print("Starting ALDEx2 to get clr-transformed values!")
    df_result_r = get_clr_function_r(df_r, cond, montecarlo_samples, instance)
    print("Finished!")
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_result = pd.DataFrame(
            ro.conversion.rpy2py(df_result_r),
            columns=metadata.index,
            index=input_counts.index,
        )
    return df_result


def get_clr(
    counts: pd.DataFrame, metadata: pd.DataFrame, mc_samples: int, r_script_path: str
) -> pd.DataFrame:
    """
    Run ALDEx2 to get CLR-transformed data. Calculates the median between all monte-carlo instances.

    Returns:
        pd.DataFrame: A dataframe with the CLR-transformed data.

    Arguments:
        counts -> pd.DataFrame: A dataframe with the raw counts.
        metadata -> pd.DataFrame: A dataframe with the metadata.
        mc_samples -> int: The number of Monte Carlo samples to use.
        r_script_path -> str: The path to the R script to run ALDEx2.
    """

    dfs = {}
    for i in range(1, mc_samples + 1):
        dfs[i] = get_clr_instance(
            counts=counts,
            metadata=metadata,
            r_script_path=r_script_path,
            instance=i,
            mc_samples=mc_samples,
        )
    averages = (
        pd.concat([each.stack() for each in dfs.values()], axis=1)
        .apply(lambda x: x.median(), axis=1)
        .unstack()
    )
    return averages


def MA_plot(
    results_input,
    effect_threshold=2,
    horizontal_line=True,
    figsize=(12, 8),
    title="MA plot",
) -> None:
    """
    Plot the MA plot of ALDEx2 results.

    Returns:
        None.

    Arguments:
        results_input -> pd.DataFrame: A dataframe with the results of ALDEx2.
        effect_threshold -> float: The threshold for the effect size value in order to mark the points red. Default is 2.
        horizontal_line -> bool: If True, a horizontal line is drawn at value 0. Default is True.
        figsize -> tuple: The size of the figure. Default is (12,8).
        title -> str: The title of the figure. Default is 'MA plot'.
    """
    results = results_input.copy()
    results["effect_hue"] = results["effect"].apply(
        lambda x: f"abs(effect) > {effect_threshold}"
        if abs(x) > effect_threshold
        else f"abs(effect) < {effect_threshold}"
    )
    results["size"] = results["effect"].apply(
        lambda x: "big" if abs(x) > effect_threshold else "small"
    )
    color_dict = dict(
        {
            f"abs(effect) > {effect_threshold}": "red",
            f"abs(effect) < {effect_threshold}": "black",
        }
    )
    size_dict = dict({"small": 5, "big": 20})
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        x="rab.all",
        y="diff.btw",
        hue="effect_hue",
        data=results,
        palette=color_dict,
        size="size",
        sizes=size_dict,
        legend="full",
        ax=ax,
    )
    if horizontal_line:
        plt.axhline(y=0, color="black", lw=1, ls="--")
    legend = ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    if len(results["size"].value_counts()) == 1:
        ax.get_legend().remove()
    else:
        ax.legend(handles=handles[1:3], labels=labels[1:3], title="")
    plt.ylabel("Difference")
    plt.xlabel("Abundance")
    plt.title(title)
    return None


def vulcano_plot(
    results_input,
    pBH_threshold=0.1,
    logFC_threshold=1,
    legend=True,
    figsize=(10, 10),
    title="Vulcano plot",
) -> None:
    """
    Plot the vulcano plot of ALDEx2 results.

    Returns:
        None.

    Arguments:
        results_input -> pd.DataFrame: A dataframe with the results of ALDEx2.
        pBH_threshold -> float: The threshold for the pBH value in order to mark the points as significant. Default is 0.1.
        logFC_threshold -> float: The threshold for the logFC value in order to mark the points as significant. Default is 1.
        legend -> bool: If True, a legend is drawn. Default is True.
        figsize -> tuple: The size of the figure. Default is (10,10).
        title -> str: The title of the figure. Default is 'Vulcano plot'.

    """
    results = results_input.copy()
    for_vulcano = (
        results[["diff.btw", "we.eBH"]].reset_index().rename(columns={"index": "gene"})
    )
    for_vulcano = for_vulcano.rename(
        columns={"diff.btw": "median_diff", "we.eBH": "pBH"}
    )

    for_vulcano["log_pValue"] = for_vulcano["pBH"].apply(lambda x: -np.log10(x))
    for_vulcano["log_FC"] = for_vulcano["median_diff"]
    logPvalues = for_vulcano.log_pValue
    logFCs = for_vulcano.log_FC
    plt.figure(figsize=figsize)

    hue = for_vulcano.gene.copy()
    hue.loc[logPvalues < -np.log10(pBH_threshold) * 0.99] = None
    ax = sns.scatterplot(y="log_pValue", x="log_FC", data=for_vulcano, hue=hue, alpha=1)
    ax = sns.scatterplot(
        y="log_pValue",
        x="log_FC",
        data=for_vulcano.loc[hue.isnull()],
        color="grey",
        marker=".",
    )

    # ax.set_xlabel('TGOLD'+'    '*10+'TOVA'+'\nlogFC', fontsize=20) # TODO: Add the labels to x axis depending on the condition names
    ax.set_xlabel("logFC", fontsize=20)
    ax_lim = np.abs(ax.get_xlim()).max()
    ax.set_xlim([-ax_lim, ax_lim])
    ax.set_ylabel("$-\log(P_{BH})$", fontsize=20)
    plt.axvline(x=-(logFC_threshold), color="red", linestyle="--", linewidth=1)
    plt.axvline(x=logFC_threshold, color="red", linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(pBH_threshold), color="red", linestyle="--", linewidth=1)
    plt.title(title)
    if legend == False:
        if ax.get_legend():
            ax.get_legend().remove()
