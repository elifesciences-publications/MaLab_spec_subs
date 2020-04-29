import pandas as pd
import numpy as np
import os

from IPython.display import display
from collections import OrderedDict
from SSanalysis import JSDcalc
from SSutility.SSerrors import SequenceAnalysisError
from SSutility.SSfasta import align_fasta_to_df

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
        test_var = col[test_species_idx][0]
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
    pre_pos = align_df.loc[idx.values[0],:align_pos-1]
    pre_pos_vc = pre_pos.value_counts()
    native_pos = align_pos
    if '-' in pre_pos_vc:
        native_pos -= pre_pos_vc['-']
    return native_pos


def calc_z_scores(scores):
    """For an array/pandas Series of scores, calculates z-scores and returns them. Ignores NaN values.

    :param scores: array/Series of scores for which z-scores will be calculated
    :return: z_scores: array/Series of same dimensions/index as scores.
    """
    mean = np.nanmean(scores)
    std = np.nanstd(scores)
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
    jsd_srs = pd.Series(index=align_df.columns)
    if keep_test_spec:
        jsd_df = align_df.copy()
    else:
        jsd_df = align_df.drop(test_idx,axis=0)
    jsd_nd = jsd_df.values
    weights = JSDcalc.calculate_sequence_weights(jsd_nd)
    for i,col in enumerate(jsd_nd.T):
        jsd = JSDcalc.JSD(col,blosum62_bg,weights,aas,use_gap_penalty=use_gap_penalty)
        # if i==42:
        #     print(use_gap_penalty)
        #     print(jsd)
        #     print("FORCE")
        #     print(col)
        #     print(JSDcalc.JSD(col, blosum62_bg, weights, aas, use_gap_penalty=False))
        #     pass
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
    """For col, return a dictionary containing information on gap_fraction,test_species/ outgroup variants, and the
    number of instances/ ratios of those variants

    :param (Series) col: Series corresponding to alignment dataframe column (including test species)
    :param test_spec_idx: index of col corresponding to test species record
    :return: metrics: dictionary containing test_variant, outgroup_variant, and gap_fraction
    """
    outgroup_variant = col.mode().iloc[0]
    gap_ratio = gap_fraction(col)
    test_variant = col[test_spec_idx][0]
    vc = col.value_counts()
    og_var_freq,test_var_freq = vc[outgroup_variant],vc[test_variant]
    og_var_ratio,test_var_ratio = og_var_freq/len(col), test_var_freq/len(col)
    vc_metrics = {'test_variant':test_variant,'outgroup_variant':outgroup_variant,'gap_fraction':gap_ratio,
               'outgroup_variant_count':og_var_freq,'test_variant_count':test_var_freq,
               'outgroup_variant_ratio':og_var_ratio,'test_variant_ratio':test_var_ratio}
    return vc_metrics

def test_outgroup_blosum(col,test_spec_idx,blos_df):
    """Calculates average BLOSUM62 score of the test species variant in col against all outgroup variants. NaN for gaps.

    :param col: pandas Series corresponding to column from align_df
    :param test_spec_idx: index corresponding to test species
    :param blos_df: DataFrame corresponding to blosum matrix; index and columns are amino acid characters
    :return: Average test variant vs outgroup variant blosum scores for col.
    """
    skip_chars = ['X','-','U']
    #test_spec_idx is Pandas Index object, use [0] to get character value instead of Series
    test_var = col[test_spec_idx][0]
    outgroup_col = col.drop(test_spec_idx)
    outgroup_col = outgroup_col[~outgroup_col.isin(skip_chars)]
    og_col_nd = outgroup_col.values
    if test_var not in skip_chars:
        blos_mean = np.mean(blos_df[test_var][og_col_nd])
    else:
        blos_mean = np.nan
    return blos_mean

def test_outgroup_blosum_series(align_df,test_spec_idx,blos_df):
    """Returns series of scores/z-scores for test vs outgroup blosum values for entire align_df.

    :param align_df: MSA DataFrame
    :param test_spec_idx: Index object corresponding to test_species. Remaining species in align_df considered outgroup
    :param blos_df: Blosum62 DataFrame
    :return: blos_srs: Series of test vs outgroup BLOSUM62 scores
    :return: blos_z: Z-scores calculated for above series
    """
    blos_srs = pd.Series(index=align_df.columns)
    for pos in align_df.columns:
        aln_col = align_df.loc[:,pos]
        blos_mean = test_outgroup_blosum(aln_col,test_spec_idx,blos_df)
        blos_srs[pos] = blos_mean
    blos_srs.name = "Test-Outgroup BLOSUM62"
    blos_z = calc_z_scores(blos_srs)
    return blos_srs, blos_z

def pairwise_outgroup_blosum(col,test_spec_idx,blos_df):
    """Returns average pairwise blosum score between all non-gap residues in outgroup of column

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

def gene_summary_table(align_df, ncbi_idx, test_idx,blos_df, display_summary=False,drop_NCBI=True,
                       summary_table_outpath="", use_jsd_gap_penalty=True):
    """Given an align_df representing OrthoDB and NCBI record multiple sequence alignment, calculates JSD, BLOSUM,
    gap and variant metrics for the OrthoDB records only.

    :param align_df: OrthoDB and NCBI MSA DataFrame. Columns are 1-indexed alignment positions and index is record ids.
    :param ncbi_idx: Record id of NCBI record. Excluded from analysis calculations if drop_NCBI=True. Stored in
    AGS Variant column in summary table
    :param test_idx: Record id of test_species (all other ODB records will be used as the outgroup for JSD and BLOSUM)
    :param blos_df: BLOSUM62 matrix as DataFrame taking amino acid chars as index/col values
    :param display_summary: If true, displays calculated summary statistics table to stdout
    :param drop_NCBI: If true, NCBI records will not be considered part of the outgroup and will be exluced from
    analysis calculations
    :return:
    """
    if drop_NCBI:
        analysis_df = align_df.drop(index=ncbi_idx)
    else:
        analysis_df = align_df
    unique_thresh = max(int(len(analysis_df)*0.1), 1)
    uniques = find_uniques(analysis_df,unique_thresh,test_idx)
    unique_pos = uniques.columns
    if len(unique_pos) == 0:
        raise SequenceAnalysisError(0,"No species unique substitutions under occurence " 
                                      "threshold {0} instances".format(unique_thresh))
    n_seq = len(analysis_df)
    #Calculate JSD and BLOSUM + z-scores for entire alignment
    jsd, jsd_z = generate_jsd_series(test_idx,analysis_df,keep_test_spec=False,use_gap_penalty=use_jsd_gap_penalty)
    test_og_blos_srs, test_og_blos_z_srs = test_outgroup_blosum_series(analysis_df,test_idx,blos_df)
    summary_col_labels = ['Test Species Position','Test Variant','AGS Variant','Test Variant Count',
                          'Outgroup Variant', 'Outgroup Variant Count', 'Analysis Sequences', 'Gap Fraction',
                          'JSD','JSD Z-Score','Test-Outgroup BLOSUM62', 'Test-Outgroup BLOSUM Z-Score',
                          'Outgroup Pairwise BLOSUM62']
    summary_df = pd.DataFrame(columns=summary_col_labels)
    summary_df.index.name = "MSA Position"
    for pos in unique_pos:
        #Calculate variant count metrics, pull blos/jsd values from alignment-calculated Series
        aln_col = uniques.loc[:,pos]
        vc_metrics = variant_counts(aln_col,test_idx)
        tv, tvc, ov, ovc, gf = [vc_metrics[label] for label in ['test_variant','test_variant_count','outgroup_variant',
                                                            'outgroup_variant_count','gap_fraction']]
        test_og_blos,test_og_blos_z = test_og_blos_srs[pos],test_og_blos_z_srs[pos]
        #Calculate outgroup pairwise BLOSUM, native sequence position of substitution in test-species
        og_pw_blos = pairwise_outgroup_blosum(aln_col,test_idx,blos_df)
        native_pos = align_pos_to_native_pos(analysis_df,test_idx,pos)
        ags_var = align_df.loc[ncbi_idx,pos][0]
        row_dict = dict(zip(summary_col_labels,[native_pos,tv,ags_var,tvc,ov,ovc,n_seq,gf,jsd[pos],
                                                jsd_z[pos],test_og_blos,test_og_blos_z,og_pw_blos]))
        summary_df.loc[pos,:] = row_dict
    if display_summary:
        print("Test Species Index: {0}".format(test_idx))
        display(summary_df)
    if summary_table_outpath:
        summary_df.to_csv(summary_table_outpath,sep='\t',float_format='%.5f')
    return summary_df

def load_summary_table(summary_fpath):
    summary_df = pd.read_csv(summary_fpath,sep='\t',index_col=0)
    return summary_df

def overall_summary_table(config, gene_symbols,use_jsd_gap_penalty=True,force_recalc=False):
    """Calculates summary analysis statistics for every gene symbol in gene_symbols.

    :param config: configparser object from config/config.txt
    :param gene_symbols: array-like or Series of gene symbols
    :param use_jsd_gap_penalty: Determines if JSD calculations are performed with gap penalty or not. If you change
    this, it is recommended to set force_recalc to True to avoid any inconsistencies between files calculated before
    the change and after.
    :param force_recalc: If True, recalculates and rewrites all summary statistic data.
    :return: None. Writes individual gene summary tables to appropriate output subdirectories and overall summary table
    to [run_name]/summary/overall_summary.tsv
    """
    from SSutility.SSerrors import load_errors, write_errors, print_errors
    run_name,errors_fname = config['RUN']['RunName'],config['RUN']['ErrorsFileName']
    ODB_test_id, NCBI_taxid = config['ODB']['ODBTestTaxID'], config['NCBI']['NCBITaxID']
    errors_fpath = "{0}/{1}".format(run_name,errors_fname)

    #Columns for overall summary table
    overall_summary_col_labels = ['Gene','MSA Position','Test Species Position', 'Test Variant', 'AGS Variant',
                                  'Test Variant Count', 'Outgroup Variant', 'Outgroup Variant Count',
                                  'Analysis Sequences', 'Gap Fraction',
                                  'JSD', 'JSD Alignment Z-Score','JSD US Z-Score',
                                  'Test-Outgroup BLOSUM62', 'Test-Outgroup BLOSUM Alignment Z-Score',
                                  'Test-Outgroup BLOSUM US Z-Score','Outgroup Pairwise BLOSUM62']
    overall_df = pd.DataFrame(columns=overall_summary_col_labels)
    overall_summary_fpath = "{0}/summary/overall_summary.tsv".format(run_name)
    check_errors, errors_df = load_errors(errors_fpath)
    for symbol in gene_symbols:
        summary_outpath = "{0}/output/{1}/{1}_summary.tsv".format(run_name,symbol)
        if not os.path.exists(summary_outpath) or force_recalc:
            #Check logged errors before attempting analysis. All logged errors will cause analysis to be skipped.
            #Logged SequenceAnalysisErrors are printed to stdout (others are passed over silently)
            if check_errors and symbol in errors_df['gene_symbol'].unique():
                sae_df = errors_df.loc[errors_df['error_type']=="SequenceAnalysisError",:]
                if symbol in sae_df['gene_symbol'].unique():
                    print_errors(sae_df,symbol)
                continue
            try:
                msa_fpath =  "{0}/output/{1}/{1}_msa.fasta".format(run_name,symbol)
                records_fpath = "{0}/output/{1}/{1}_records.tsv".format(run_name,symbol)
                records_df = pd.read_csv(records_fpath,sep='\t',index_col=0)
                records_df.index.name = "record_id"
                ncbi_idx = records_df.loc[records_df['db_source']=="NCBI",:].index
                test_idx = records_df.loc[records_df['organism_taxid']==ODB_test_id,:].index
                align_df = align_fasta_to_df(msa_fpath)

                summary_df = gene_summary_table(align_df,ncbi_idx,test_idx,blos_df,
                               display_summary=False,drop_NCBI=True,summary_table_outpath=summary_outpath,
                                use_jsd_gap_penalty=use_jsd_gap_penalty)
            except SequenceAnalysisError as sae:
                write_errors(errors_fpath,symbol,sae)
                continue
        else:
            summary_df = load_summary_table(summary_outpath)
        #Format summary_df into overall_summary format (add Gene and MSA position columns, rename Z-score columns)
        formatted = summary_df.reset_index(drop=False)
        formatted.insert(0,"Gene",[symbol]*len(formatted))
        formatted = formatted.rename(columns={'JSD Z-Score':'JSD Alignment Z-Score',
                                       'Test-Outgroup BLOSUM Z-Score':'Test-Outgroup BLOSUM Alignment Z-Score'})
        overall_df = overall_df.append(formatted,ignore_index=True,sort=False)
    display_overall=False
    if display_overall:
        with pd.option_context('display.max_columns',None):
            display(overall_df)
    #Unique Substitution Set-wide Z scores for JSD and BLOSUM
    us_jsd, us_blos = calc_z_scores(overall_df['JSD']), calc_z_scores(overall_df['Test-Outgroup BLOSUM62'])
    overall_df.loc[:,'JSD US Z-Score'] = us_jsd
    overall_df.loc[:,'Test-Outgroup BLOSUM US Z-Score'] = us_blos
    overall_df.to_csv(overall_summary_fpath,sep='\t',float_format='%.5f')
    return overall_df


#Global variables for module
aas, blosum62_bg, blos_df, sim_matrix = gen_blos_df()