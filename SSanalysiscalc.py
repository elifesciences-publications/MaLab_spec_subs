import pandas as pd
from IPython.display import display
from collections import OrderedDict
from Bio.SubsMat.MatrixInfo import blosum62

def gen_blos_df():
    """Background distribution data from published code from below citation
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    :return:
    """
    global aas, blosum62_bg
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

def find_uniques(align_df, sub_freq_threshold, test_species_idx,display_uniques):
    #align
    uniques = pd.DataFrame(index=align_df.index)
    for pos in align_df:
        col = align_df[pos]
        sub_counts = col.value_counts()
        c_counts = sub_counts.value_counts()
        for index,value in sub_counts.iteritems():
    #         if(value == 1):# and c_counts.loc[value] == 1:
            if value <= sub_freq_threshold:
                if col[test_species_idx] == index and index != 'X':
                    uniques.loc[:,pos] = col
                    break
    if display_uniques:
        print("Threshold Number of Sequences: "+str(int(sub_freq_threshold)))#/len(ordered)))
        display(uniques)
    return uniques

def gap_fraction(col):
    return sum([c == '-' for c in col]) / len(col)

def sub_fraction(col, sub):
    return sum([c == sub for c in col]) / len(col)

def species_gaps_dict(align_df, gap_char='-'):
    gaps_dict = OrderedDict()
    unique_gaps = OrderedDict()
    unique_subs = OrderedDict()
    align_nd = align_df.values
    for j, spec in enumerate(align_df.index):
        gaps_dict[spec] = []
        unique_gaps[spec] = []
        unique_subs[spec] = []
        spec_row = align_df.loc[spec, :].values
        #         row = align_nd[j]
        for i, entry in enumerate(spec_row):
            #             outgroup_col_ = align_df.iloc[:,i].copy()
            outgroup_col_ = align_nd.T[i].copy()
            #             outgroup_col_.drop(spec,inplace=True)
            outgroup_col_ = np.delete(outgroup_col_, obj=j)
            ratio = gap_fraction(outgroup_col_)
            if entry == gap_char:
                gaps_dict[spec].append(i + 1)
                # Check if remainder of positions have gap; if none do, add to unique gaps.
                if ratio <= 0.1:
                    unique_gaps[spec].append(i + 1)
            else:
                if ratio >= 0.9:
                    unique_subs[spec].append(i + 1)

    return gaps_dict, unique_gaps, unique_subs



