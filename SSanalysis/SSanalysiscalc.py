import pandas as pd
import numpy as np
from IPython.display import display
from collections import OrderedDict
from SSanalysis import JSDcalc

def gen_blos_df():
    from Bio.SubsMat.MatrixInfo import blosum62
    """Background distribution data from published code from below citation
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    :return: aas: list of amino acid chars inlcuding gap char '-'
    :return: blosum62_bg: background distribution of amino acids from BLOSUM dataset
    :return: blos_df: DataFrame corresponding to BLOSUM62 matrix, index/cols are amino acid characters
    :return sim_matrix: np.ndarray of values in blos_df
    """
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-']
    bg_probs = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, \
                0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
    blosum62_bg = bg_probs
    blos_df = pd.DataFrame(index=aas[:-1],columns=aas[:-1])
    for pair in blosum62:
        val = blosum62[pair]
        first, second = pair[0],pair[1]
        if first in aas and second in aas:
            blos_df.loc[first,second] = val
            blos_df.loc[second,first] = val
    sim_matrix = blos_df.values
    return aas, blosum62_bg, blos_df, sim_matrix

def find_uniques(align_df, sub_freq_threshold, test_species_idx,display_uniques=False):
    """Identifies unnique positions from align_df with frequency <= sub_fre_threshold

    :param align_df: DataFrame of alignment characters. Columns: 1-indexed alignment positions, Index: record ids
    :param (int) sub_freq_threshold: max allowed number of instances of a substitution for it to be considered unique
    :param test_species_idx: record_id for test_species record in align_df for which uniques will be identified
    :param (boolean) display_uniques: If true, displays table of unique residues identified
    :return:
    """
    uniques = pd.DataFrame(index=align_df.index)
    for pos in align_df:
        col = align_df[pos]
        sub_counts = col.value_counts()
        test_var = col[test_species_idx]
        test_vc = sub_counts[test_var]
        if test_vc <= sub_freq_threshold and test_var!='X':
            uniques.loc[:,pos] = col
    if display_uniques:
        print("Threshold Number of Sequences: "+str(int(sub_freq_threshold)))#/len(ordered)))
        display(uniques)
    return uniques

def gap_fraction(col):
    #Fraction of positions in col corresponding to gap character '-'
    return sum([c == '-' for c in col]) / len(col)

def sub_fraction(col, sub):
    #Fraction of positions in col corresponding to character sub.
    return sum([c == sub for c in col]) / len(col)

def align_pos_to_native_pos(align_df, idx, align_pos):
    """Converts alignment position to native sequence position (ie ignores gaps).

    :param align_df: DataFrame of msa characters
    :param idx: Index (row) of align_df to use for native position
    :param align_pos: alignment column position to convert to native position
    :return: native_pos: converted alignment position
    """
    pre_pos = align_df.loc[idx,:align_pos-1]
    pre_pos_vc = pre_pos.value_counts()
    native_pos = align_pos
    if '-' in pre_pos_vc:
        native_pos -= pre_pos_vc['-']
    return native_pos


def calc_z_scores(scores):
    """For an array/pandas Series of scores, calculates z-scores and returns them.

    :param scores: array/Series of scores for which z-scores will be calculated
    :return: z_scores: array/Series of same dimensions/index as scores.
    """
    mean = np.mean(scores)
    std = np.std(scores)
    z_scores = (scores-mean)/std
    return z_scores


def generate_jsd_series(test_idx,align_df,
                        keep_test_spec=False,use_gap_penalty=True):
    """ Generate a Jensen Shannon Difference Series for position indices within the MSA. JSD is by default based on
    align_df with test_species dropped (ie outgroup species). See JSDcalc.py for calculation details. Calculations
    based on below citation.

    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param test_species: index value in align_df corresponding to test_species. Row will be removed from calculation
    unless keep_test_spec set to True.
    :param align_df: DataFrame containing multiple sequence alignment (1-indexed positions as columns, record ids as index
    :param (boolean) keep_test_spec: Determines whether JSD calculation is done on outgroup or all records in align_df
    :return jsd_srs: pandas Series containing JSD values at each position (1-indexed)
    :return jsd_zscores: pandas Series containing JSD z-score (calculated using mean and std of this alignment only)
    """
    from SSanalysis import aas, blosum62_bg
    jsd_srs = pd.Series(index=align_df.columns)#range(1,alignlen))
    if keep_test_spec:
        jsd_df = align_df.copy()
    else:
        jsd_df = align_df.drop(test_idx,axis=0)
    jsd_nd = jsd_df.values
    weights = JSDcalc.calculate_sequence_weights(jsd_nd)
    for i,col in enumerate(jsd_nd.T):
        jsd = JSDcalc.JSD(col,blosum62_bg,weights,aas,use_gap_penalty=use_gap_penalty)
        jsd_srs[i+1] = jsd
    jsd_zscores = calc_z_scores(jsd_srs)
    return jsd_srs, jsd_zscores

def filter_score_srs(score_srs, uniques):
    """Filters score series down to positions in uniques.

    :param score_srs: pandas Series of score entries, indexed on alignment positions
    :param uniques: DataFrame corresponding to columns from alignment DataFrame with unique residue substitutions
    :return: unique_scores: pandas Series corresponding to entries in score_srs corresponding to alignment positions
    in uniques
    """
    unique_positions = uniques.columns
    unique_scores = score_srs.loc[unique_positions]
    return unique_scores

def variant_counts(col,test_spec_idx):
    """

    :param col: Alignment DataFrame column (including test species)
    :param test_spec_idx: index of col corresponding to test species record
    :return: metrics: dictionary containing test_variant, outgroup_variant, and gap_fraction
    """
    align_col = col.values
    outgroup_variant = align_col.mode().iloc[0]
    gap_ratio = JSDcalc.gap_fraction(align_col)
    test_variant = col[test_spec_idx]#align_df.loc[test_species, pos]
    vc = col.value_counts()
    og_var_freq = vc[outgroup_variant]
    test_var_freq = vc[test_variant]
    metrics = {'test_variant':test_variant,'outgroup_variant':outgroup_variant,'gap_fraction':gap_ratio,
               'outgroup_variant_freq':og_var_freq,'test_variant_freq':test_var_freq}
    return metrics

def test_outgroup_blosum(col,test_spec_idx,blos_df):
    """Calculates average BLOSUM62 score of the test species variant in col against all outgroup variants. NaN for gaps.

    :param col: pandas Series corresponding to column from align_df
    :param test_spec_idx: index corresponding to test species
    :param blos_df: DataFrame corresponding to blosum matrix; index and columns are amino acid characters
    :return: Average test variant vs outgroup variant blosum scores for col.
    """
    test_var = col[test_spec_idx]
    outgroup_col = col.drop(test_spec_idx)
    og_col_nd = outgroup_col.values
    if test_var != '-':
        blos_mean = np.mean(blos_df[test_var][og_col_nd])
    else:
        blos_mean = np.nan
    return blos_mean

def pairwise_outgroup_blosum(col,test_spec_idx,blos_df):
    """Returns average pairwise blosum score between all residues in outgroup of column

    :param col: pandas Series corresponding to column from alignment DataFrame
    :param test_spec_idx: index of test_species in col
    :param blos_df: BLOSUM DataFrame from gen_blos_df
    :return: Average pairwise BLOSUM score over all residues in outgroup against each other
    """
    outgroup_col = col.drop(test_spec_idx)
    og_col = outgroup_col.values
    pairwise_scores = []
    skip_chars = ['-', 'X']
    for i, first in enumerate(og_col):
        for j, second in enumerate(og_col):
            if i < j and not (first in skip_chars) and not second in skip_chars:
                pairwise_scores.append(blos_df[first][second])
            else:
                pass
    if pairwise_scores:
        blos_pw = np.mean(pairwise_scores)
    else:
        blos_pw = np.nan
    return blos_pw

def summary_table(align_df, uniques, test_idx,blos_df, display_summary):

    jsd, jsd_z = generate_jsd_series(test_idx,align_df)
    jsd_uniques, jsd_z_uniques = filter_score_srs(jsd,uniques), filter_score_srs(jsd_z,uniques)

    jsd_nd = align_df.copy().drop(index=test_species).values
    unique_positions = sw_metrics.index
    unique_positions.name = "MSA Position"
    n_seq = len(align_df.index)
    cols_list = ["Test Species Position", "Test Variant","Test Variant Instances","Outgroup Variant",\
                 "Outgroup Variant Instances","Aligned Sequences", "Gap Fraction", "JSD","JSD Z-Score", \
                 "Test vs Outgroup Blosum62", "Outgroup Pairwise Blosum62"]
    summary_df = pd.DataFrame(index=unique_positions,columns=cols_list)
    for pos in unique_positions:
        native_pos = JSDcalc.align_pos_to_native_pos(align_df, test_species,pos)
        jsd = pos_jsds[pos]
        jsd_z = pos_jsdzs[pos]
        align_col_srs = align_df.loc[:,pos]
        align_col = align_col_srs.values
        outgroup_variant = align_col_srs.mode().iloc[0]
        gap_fraction = JSDcalc.gap_ratio(align_col)
        test_variant = align_df.loc[test_species,pos]
        vc = align_df[pos].value_counts()
        og_var_freq = vc[outgroup_variant]
        test_var_freq = vc[test_variant]
        col = jsd_nd.T[pos-1]
        if not test_variant == '-':
            blos_mean = np.mean(blos_df[test_variant][col])
        else:
            blos_mean = np.nan
        #pairwise calculations
        others_col = jsd_nd.T[pos-1]
        pairwise_scores = []
        skip_chars = ['-','X']
        for i,first in enumerate(others_col):
            for j,second in enumerate(others_col):
                if i<j:
                    if not first in skip_chars and not second in skip_chars:
                        pairwise_scores.append(blos_df[first][second])
                    else:
                        pass
        if pairwise_scores:
#             blosum62_pairwise[pos] = np.mean(pairwise_scores)
            blos_pw = np.mean(pairwise_scores)
        else:
            blos_pw = np.nan
        row_values = [native_pos, test_variant, test_var_freq,outgroup_variant,og_var_freq,n_seq,gap_fraction,jsd,jsd_z,blos_mean,blos_pw]
        summary_row = pd.Series(name=pos,data=dict(zip(cols_list,row_values)))
#         summary_df = summary_df.append(summary_row)
        summary_df.loc[pos] = summary_row
    #Output conditional on config option
    if display_summary:
        print(test_species+" Semi-Unique Substitutions Summary")
        print("Number of Species: "+str(len(align_df.index)))
        with pd.option_context("display.max_rows",None,"display.max_columns", None,"display.max_colwidth",200):
            display(summary_df)
#             pass
    return summary_df