import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib as mpl

from IPython.display import display

# mpl.rcParams['figure.dpi'] = 200
# mpl.rcParams['savefig.dpi'] = 200

def read_overall_summary(config):
    overall_fpath = "{0}/summary/overall_summary.tsv".format(config['RUN']['RunName'])
    overall_df = pd.read_csv(overall_fpath,sep='\t',index_col=0)
    return overall_df

def get_scatter_alpha(overall_df):
    #Returns alpha value to use for scatter plots using overall_df
    alpha = 385 / len(overall_df)
    if alpha < 0.1:
        alpha *= 1.3
    return alpha

def AGS_13LGS_consensus(overall_df):
    """Returns DataFrame for which AGS Variant and Test Variant are equal. Assumes Test Variant is

    :param overall_df: DataFrame of run_wide unique substitutions
    :return: gs_consensus_df: DataFrame of subset of rows from overall_df where
    """
    gs_consensus_df = overall_df.loc[overall_df['Test Variant']==overall_df['AGS Variant'],:]
    return gs_consensus_df


def gene_split(genename, summary_df):
    gene_df = summary_df[summary_df["Gene"] == genename]
    rest_df = summary_df[~ (summary_df["Gene"] == genename)]
    return gene_df, rest_df


def gene_split_scatter(genename, out_dir, title_str, summary_df, use_zscore=False):
    plt.figure()

    gene_df, rest_df = gene_split(genename, summary_df)
    if use_zscore:
        x_field = "JSD Z-Score"
        #         title = "Distribution of JSD Z-Score/ Blosum for all input unique subs", "+genename+" in red""
        title = "Distribution of JSD Z-Score/ BLOSUM, all cDNA screen hits"
        figure_path = out_dir + genename + "_Z-score_scatter.png"
    else:
        x_field = "JSD"
        #         title = "Distribution of JSD/ Blosum for all input unique subs, "+genename+" in red"
        figure_path = out_dir + genename + "_rawJSD_scatter.png"
    plt.scatter(x=rest_df[x_field], y=rest_df["Test vs Outgroup Blosum62"], alpha=alpha)
    plt.scatter(x=gene_df[x_field], y=gene_df["Test vs Outgroup Blosum62"], c="red", alpha=1)
    plt.xlabel("Position " + x_field)
    #     plt.ylabel("Average Test Species vs Outgroup Blosum62")
    plt.ylabel(SCATTER_YLABEL)
    plt.title(title_str)
    plt.savefig(figure_path)

def main():
    os.chdir("..")
    from SSutility import config
    overall_df = read_overall_summary(config)
    gs_consensus = AGS_13LGS_consensus(overall_df)
    # display(gs_consensus)
    count = 0
    for i,row in overall_df.iterrows():#gs_consensus.iterrows():
        if row['Test Variant'] == row['AGS Variant']:
            count += 1
    print(count)

if __name__ == "__main__":
    main()