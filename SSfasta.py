from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import pandas as pd

# Fasta file reading functions:
# filter_fasta_infile reads input files and outputs all records corresponding to filtered_ids to a new file
# Remaining functions provide conversions between fasta files, pandas Series, and pandas dataframes
# having alignment positions as columns

def filter_fasta_infile(filtered_ids, infile_path, outfile_path=None, ordered=False):
    # If outfile_path is provided, write filtered fasta to outfile_path
    """Generates new fasta file to outfile_path using the subset of sequences in infile_path
    which have ids in filtered_ids
    ordered: if true, sequences will be returned/ written in order of filtered_ids
             if false, uses sequence order of sequences in infile_path
    """

    def filtered_generator(filtered_ids, infile_path):
        fasta_seqs = SeqIO.parse(open(infile_path), 'fasta')
        for fasta in fasta_seqs:
            if fasta.id in filtered_ids:
                yield fasta

    def ordered_filtered_generator(filtered_ids, infile_path):
        for id_ in filtered_ids:
            fasta_seqs = SeqIO.parse(open(infile_path), 'fasta')
            for fasta in fasta_seqs:
                if fasta.id == id_:
                    yield fasta
                    break

    if outfile_path:
        if ordered:
            filtered = ordered_filtered_generator(filtered_ids, infile_path)
        else:
            filtered = filtered_generator(filtered_ids, infile_path)
        SeqIO.write(filtered, outfile_path, "fasta")
    if ordered:
        filtered = ordered_filtered_generator(filtered_ids, infile_path)
    else:
        filtered = filtered_generator(filtered_ids, infile_path)
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
    fasta_seqs = SeqIO.parse(open(fasta_path), 'fasta')
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
    dist_srs = avg_dist_srs(align_srs.index, id_dm)
    return align_df, dist_srs


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
