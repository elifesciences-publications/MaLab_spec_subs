from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import pandas as pd
import numpy as np
import subprocess
import warnings
import os
import SSerrors
# Fasta file reading functions:
# filter_fasta_infile reads input files and outputs all records corresponding to filtered_ids to a new file
# Remaining functions provide conversions between fasta files, pandas Series, and pandas dataframes
# having alignment positions as columns

def _filtered_generator(filtered_ids, infile_path):
    """Generator function, yields Bio.Seq fastas from infile_path in order of appearance in infile"""
    fasta_f = open(infile_path)
    fasta_seqs = SeqIO.parse(fasta_f, 'fasta')
    for fasta in fasta_seqs:
        if fasta.id in filtered_ids:
            yield fasta
    fasta_f.close()

def _ordered_filtered_generator(filtered_ids, infile_path):
    """Generator function, yields Bio.Seq fastas from infile_path in order of appearance in filtered_ids"""
    fasta_f = open(infile_path)
    for id_ in filtered_ids:
        fasta_seqs = SeqIO.parse(fasta_f, 'fasta')
        for fasta in fasta_seqs:
            if fasta.id == id_:
                yield fasta
                break
    fasta_f.close()

def _generator_wrapper(filtered_ids, infile_path, ordered, missing_warning=False):
    """Wrapper function which generates appropriate generatoe object based on ordered. Issues warnings for missing
    values from filtered_ids if not present in infile_path if missing_warning is True.

    :param filtered_ids: See filter_fasta_infile
    :param infile_path: See filter_fasta_infile
    :param ordered: if True, generator returns records in order of appearance in filtered_ids; False -> order from
    record order in infile.
    :param missing_warning: boolean on whether to issue warnings for missing filtered_id values from infile
    :return:
    """
    if missing_warning:
        with open(infile_path) as infile_f:
            infile_ids = [fasta.id for fasta in SeqIO.parse(infile_f, 'fasta')]
            for record_id in filtered_ids:
                if record_id not in infile_ids:
                    msg = "Infile {0} is missing record id {1} from filtered_ids".format(infile_path, record_id)
                    warnings.warn(msg)
    if ordered:
        filtered = _ordered_filtered_generator(filtered_ids, infile_path)
    else:
        filtered = _filtered_generator(filtered_ids, infile_path)
    return filtered
def filter_fasta_infile(filtered_ids, infile_path, outfile_path=None, ordered=False):
    """Filters fasta_infile records from infile_path and returns a series of fasta records if id in filtered_ids

    :param filtered_ids: Iterable containing all id values for which records will be saved
    :param infile_path: Fasta file containing unfiltered set of records. If missing records from filtered_ids,
    will issue warning
    :param outfile_path: Optional filepath parameter. If provided, filtered results will be written to file path
    :param ordered: boolean, if true, sequences will be returned/ written in order of filtered_ids
             if false, uses sequence order of sequences in infile_path
    :return: Series of fasta sequences (index is fasta.id) for which id is present in filtered_ids
    """
    if outfile_path:
        filtered = _generator_wrapper(filtered_ids,infile_path,ordered)
        SeqIO.write(filtered, outfile_path, "fasta")
    with warnings.catch_warnings(record=True) as w:
        filtered = _generator_wrapper(filtered_ids, infile_path, ordered, True)
    filtered_srs = pd.Series(index=filtered_ids)
    for fasta in filtered:
        filtered_srs[fasta.id] = str(fasta.seq)
    return filtered_srs


def srs_to_fasta(seq_srs, outfile_path):
    # Write records in seq_srs to outfile_path in fasta format
    def record_generator(seq_srs):
        for idx, seq in seq_srs.iteritems():
            record = SeqRecord(Seq(seq, IUPAC.protein), id=idx)
            yield record

    records = record_generator(seq_srs)
    SeqIO.write(records, outfile_path, "fasta")


def fasta_to_srs(fasta_path):
    with open(fasta_path) as fasta_f:
        fasta_seqs = SeqIO.parse(fasta_f, 'fasta')
        id_seq_map = OrderedDict()
        for fasta in fasta_seqs:
            record_id = fasta.id
            seq = str(fasta.seq)
            id_seq_map[record_id] = seq
        return pd.Series(name="seq", data=id_seq_map)


def align_srs_to_df(align_srs):
    # Returns DataFrame object from series of aligned sequences; columns are 1-indexed positions
    # Values are characters in alignment, index is ODB sequence IDs
    n_seq = len(align_srs)
    #     display(align_srs)
    #     display(align_srs.iloc[0])
    seq_len = len(align_srs.iloc[0])
    align_df = pd.DataFrame(index=align_srs.index, columns=range(seq_len))
    for idx, seq in align_srs.iteritems():
        align_df.loc[idx, :] = list(seq)
    align_df.columns += 1
    return align_df


def seq_srs_to_align_df(seq_srs, align_in_fpath, align_out_fpath):
    """Transform seq_srs (pandas Series containing sequence texts) to a DataFrame for which each column
    is an alignment position and column. Writes input fasta and output fastas for alignment to align_in_fpath
    and align_out_fpath respectively. Also returns average (non-diagonal) identity distances"""
    srs_to_fasta(seq_srs, align_in_fpath)
    n, ordered_ids, id_dm, align_srs = construct_id_dm(seq_srs, align_in_fpath, align_out_fpath)
    align_df = align_srs_to_df(align_srs)
    # dist_srs = avg_dist_srs(align_srs.index, id_dm)
    return align_df#, dist_srs


def align_srs_to_seq_srs(align_srs, outfile_path=None):
    # Return new Series (same index) of sequences with gap characters dropped
    # If outfile_path is provided, write un-aligned record seqs to new fasta file
    seq_srs = pd.Series(index=align_srs.index)
    for idx, align_seq in align_srs.iteritems():
        seq = align_seq.replace("-", "")
        seq_srs[idx] = seq
    if outfile_path:
        srs_to_fasta(seq_srs, outfile_path)
    return seq


def align_df_to_srs(align_df):
    # Returns series of aligned sequences from array of aligned positions
    align_srs = pd.Series(index=align_df.index)
    for idx, record in align_df.iterrows():
        #       #seq is a string joining all characters with no delimiter (i.e. the original aligned sequence with gaps)
        seq = ''.join(record.values)
        align_srs[idx] = seq
    return align_srs

### Record DataFrame Functions
def load_tsv_table(input_tsv_fpath,tax_subset=[],ODB_ID_index=True):
    """Loads OrthoDB tsv data into a pandas DataFrame

    :param input_tsv_fpath: File path to OrthoDB tsv.
    :param (array-like) tax_subset: If provided, returned DataFrame only contains records which have organism_taxid
    in tax_subset
    :param (boolean) ODB_ID_index: If True, index on int_prot_id (OrthoDB record identifier string), else int-indexed
    :return: tsv_df: DataFrame containing records from input_tsv_fpath with above filters.
    """
    if not os.path.exists(input_tsv_fpath):
        raise SSerrors.RecordDataError(0,"Missing File at path: {0}".format(input_tsv_fpath))
    tsv_df = pd.read_csv(input_tsv_fpath, delimiter='\t')
    if ODB_ID_index:
        tsv_df = tsv_df.set_index(keys="int_prot_id", drop=True)  # drop=False)
    if len(tax_subset) > 0:
        tsv_df = tsv_df.loc[tsv_df["organism_taxid"].isin(tax_subset), :]
    return tsv_df


    ### Distance Matrix Functions
#construct_id_dm makes an np.ndarray for the identity distance matrix of sequences for which OrthoDB id is
# in the index of seq_df; distance matrix rows will be ordered to the order in seq_fpath if ordered isFalse
def construct_id_dm(seq_df, seq_fpath, align_outpath="tmp/iddm_align.fasta", ordered=False,
                    kalign_silent=True):
    """Constructs an np.ndarray corresponding to the identity distance matrix of records in seq_df

    :param seq_df: DataFrame of OrthoDB/ NCBI sequence records; should only contain records for which identity
    distance matrix will be computed
    :param seq_fpath:  Path of fasta file containing at least all of the records in seq_df. Can contain more records
    than are in seq_df - a temporary file containing only the records in seq_df.index will be generated (filtered_fpath)
    :param align_outpath: Optional filepath. If provided, the resulting alignment will be stored there. Otherwise,
    written to a temporary file (tmp/iddm_align.fasta)
    :param ordered: boolean. True: distance matrix rows will be ordered by the order of records in seq_df.index;
    False: distance matrix rows will be ordered by the order of records in seq_fpath
    :return: id_dm: np.ndarray of identity distance matrix calculated by AlignIO
    :return: align_srs: pandas Series object containing aligned sequences
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio import AlignIO
    # Filter records in seq_fpath to new fasta only containing records in seq_df.index
    # filtered_outpath = "tmp/iddm.fasta"
    filtered_fpath = "tmp/alias_matches.fasta"
    filter_fasta_infile(seq_df.index, seq_fpath, outfile_path=filtered_fpath, ordered=ordered)
    # KAlign sequences in filtered_outpath, write to align_outpath
    # n, ordered_ids, ka_dm, align_outfile = load_ka_distmat(filtered_outpath, align_outfile=align_outpath)
    # proc = subprocess.run(args=["kalign", '-i', filtered_fpath, "-o", align_outpath, "-f", "fasta"])
    #Use subprocess.Popen instead of subprocess.run for PyCharm compatability
    from subprocess import Popen
    with open(filtered_fpath,'r') as filtered_f, open(align_outpath,'wt',encoding='utf-8') as align_f:
        args = ['kalign']
        if kalign_silent:
            subprocess.run(args=args, stdin=filtered_f, stdout=align_f, stderr=subprocess.PIPE, text=True)
        else:
            subprocess.run(args=args, stdin=filtered_f, stdout=align_f, text=True)

    align_srs = fasta_to_srs(align_outpath)
    with open(align_outpath) as aligned_f:
        aln = AlignIO.read(aligned_f, 'fasta')
    calculator = DistanceCalculator('identity')
    id_dm_obj = calculator.get_distance(aln)
    # Convert AlignIO object to np.ndarray
    for i, r in enumerate(id_dm_obj):
        if i == 0:
            id_dm = np.array(r)
        else:
            id_dm = np.vstack((id_dm, r))
    return id_dm, align_srs

def avg_dist_srs(index,distmat):
    #index is a pandas Index object with entries corresponding to the distmat (i.e. lengths should be equal)
    #Calculate mean of non-self record distances (diagonal distances generally force-set to 0, so
    #sum functions as intended)
    n = len(distmat)
    avg_dists = np.sum(distmat, axis=1)/(n-1)
    dist_srs = pd.Series(data=avg_dists,index=index,name="dist")
    return dist_srs