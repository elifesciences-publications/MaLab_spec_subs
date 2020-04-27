import pandas as pd
import numpy as np
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


# Calculates Henikoff (1994) sequence weights form msa
from collections import Counter


def calculate_sequence_weights(msa, aas):
    """ Calculate the sequence weights using the Henikoff '94 method
    for the given msa. """
    n = len(msa)
    alignlen = len(msa[0])
    seq_weights = np.zeros((n))
    for col in msa.T:
        freqs = Counter(col)
        del freqs['-']
        num_observed = len(freqs)

        #         d_ = [freqs[aa]*num_observed if aa != '-' else 0 for aa in col]
        d_ = [freqs[aa] * num_observed for aa in col]
        inverse = [1 / d if d > 0 else 0 for d in d_]
        seq_weights += inverse
    seq_weights /= alignlen
    return seq_weights


def weighted_fcpseudo(col, seq_weights, pc_amount=0.000001):
    """Generates Pc distribution for column col. Distribution is weighted using seq_weights

    :param col:
    :param seq_weights:
    :param pc_amount:
    :return:
    """
    if len(seq_weights) != len(col):
        return [1]*len(col)
    else:
        freqs = np.array([pc_amount]*len(aas))
        for i,aa in enumerate(aas):
#             print(col)
#             print([seq_weights[j] if aa == entry else 0 for j,entry in enumerate(col)])
            freqs[i] += sum([seq_weights[j] if aa == entry else 0 for j,entry in enumerate(col)])
        for i,entry in enumerate(freqs):
#             if entry > pc_amount:
#                  freqs[i] -= pc_amount
            freqs[i] /= (sum(seq_weights) + len(aas)*pc_amount)
        return freqs


def weighted_gap_penalty(col,weights):
    #Return simple gap penalty multiplier for column
    gap_weights = np.array([w if c == "-" else 0 for c,w in zip(col,weights)])
    gap_sum = sum(gap_weights)
    return 1 - (gap_sum/sum(weights))
# print(weighted_gap_penalty(['A','A','-'],[0.3333,0.3333,0.3333]))
def relative_entropy(pC,bg,gap_penalty=0):
    """Calculate the relative entropy of the column distribution with a
    background distribution specified in bg_distr. This is similar to the
    approach proposed in Wang and Samudrala 06."""
    if len(bg) == 20 and len(pC) == 21:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    ratio = pC/bg
    log_ratio = np.log(ratio)/np.log(2)
    RE = sum([a*b for a,b in zip(pC,log_ratio)])
    if gap_penalty:
        RE*= gap_penalty
    return RE

def col_RE(col,bg,weights,gap_penalty=1):
    pC = weighted_fcpseudo(col,weights)
    if len(bg) == 20 and len(pC) == 21:
        pC = pC[:-1]
        pC /= sum(pC)
    ratio = pC/bg
    log_ratio = np.log(ratio)/np.log(2)
    RE = sum([a*b for a,b in zip(pC,log_ratio)])
    if gap_penalty:
        RE*= gap_penalty
    return RE
# relative_entropy(np.array([0.25,0.25,0.25,0.25]),np.array([0.1,0.1,0.1,0.7]))

def JSD(col,bg,weights,gap_penalty=1,jsd_lam=0.5):
    pC = weighted_fcpseudo(col,weights)
    if len(bg) == 20:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    gap_penalty = weighted_gap_penalty(col,weights)
    r_dist = np.array([jsd_lam*aa for aa in pC]) + np.array([(1-jsd_lam)*aa for aa in bg])

    RE_pCr = relative_entropy(pC,r_dist,gap_penalty)
    RE_qr = relative_entropy(bg,r_dist, gap_penalty)
    jsd = jsd_lam*RE_pCr + (1-jsd_lam)*RE_qr
    return jsd

def calc_z_scores(scores):
    """

    :param scores:
    :return:
    """
    mean = np.mean(scores)
    std = np.std(scores)
    z_scores = (scores-mean)/std
    return z_scores

#Generate a Jensen Shannon Difference Series for position indices within the MSA
def generate_jsd_series(test_species,align_df, alignlen, aas,blosum62_bg):
    jsd_srs = pd.Series(index=range(1,alignlen))
    jsd_lam = 0.5
    gap_threshold = 0.3
    gap_positions = []
    jsd_df = align_df.drop(test_species,axis=0)
    # jsd_df = align_df.copy()
    jsd_nd = jsd_df.values
    weights = calculate_sequence_weights(jsd_nd,aas)
    # print(weights)
    for i,col in enumerate(jsd_nd.T):
        gap_penalty = weighted_gap_penalty(col,weights)
        jsd = JSD(col,blosum62_bg,weights,gap_penalty)
        jsd_srs[i+1] = jsd
        if gap_fraction(col) > gap_threshold:
            gap_positions.append(i+1)
    jsd_zscores = calc_z_scores(jsd_srs)
    with pd.option_context("display.max_rows",None,"display.max_columns", None):
#         display(jsd_srs)
        pass
    return jsd_srs, jsd_zscores, gap_positions