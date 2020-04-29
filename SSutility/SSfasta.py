#SSfasta.py - Fasta file managing/ conversion/ filtering functions
#and logging helper functions for storing errors.
# Copyright (C) 2020  Evan Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


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
from SSutility import SSerrors

###Record filtering functions###

def ordered_record_generator(fpath,ordered_ids):
    """Generates Seq records from fpath, limited to ordered_ids if provided. Generator order is order in ordered_ids.
    Internally tracks yielded record ids to eliminate duplicates.
    :param fpath: Fasta file path
    :param ordered_ids: array-like containing ordered record ids
    :return:
    """
    yielded = set()
    for id in ordered_ids:
        with open(fpath) as fasta_f:
            fastas = SeqIO.parse(fasta_f,'fasta')
            for fasta in fastas:
                if fasta.id == id and fasta.id not in yielded:
                    yielded.add(fasta.id)
                    yield fasta
                    break

def record_generator(fpath,ids=[]):
    """Generates Seq records from fpath, limited to ids if provided. Ordered as in fpath
    :param fpath: Fasta file path
    :param ids: If provided, only yield records corresponding to record ids provided. Default value causes all records
    to be present in generator
    :return: generator object
    """
    yielded = set()
    with open(fpath) as fasta_f:
        fastas = SeqIO.parse(fasta_f, 'fasta')
        for fasta in fastas:
            if ((len(ids)>0 and fasta.id in ids) or \
                    (len(ids)==0)) and fasta.id not in yielded:
                yielded.add(fasta.id)
                yield fasta

def ODB_NCBI_generator(ODB_fpath,NCBI_fpath,odb_subset=[],ncbi_subset=[],ordered=False):
    """Generator object that yields all Bio.Seq objects in ODB_fpath and then NCBI_fpath, filtered down with optional
    parameters odb_subset and ncbi_subset.

    :param ODB_fpath: file path to OrthoDB fasta
    :param NCBI_fpath: file path to NCBI fasta
    :param odb_subset: If provided, only ODB records present in odb_subset will be provided by generator. Else all
    records will be provided
    :param ncbi_subset: If provided, only NCBI records present in ncbi_subset will be provided by generator. Else all
    records will be provided
    :param ordered: If true, odb_subset and ncbi_subset must be provided. Causes records to be returned in order they
    are present in odb_subset and ncbi_subset.
    :return: Generator object which provides Bio.Seq objects filtered and ordered as described above.
    """

    if ordered:
        if len(odb_subset) == 0 or len(ncbi_subset)==0:
            raise ValueError("Both odb_subset and ncbi_subset must be provided if ordered is True.")
        else:
            odb_generator = ordered_record_generator(ODB_fpath,odb_subset)
            ncbi_generator = ordered_record_generator(NCBI_fpath,ncbi_subset)
    else:
        odb_generator = record_generator(ODB_fpath, odb_subset)
        ncbi_generator = record_generator(NCBI_fpath, ncbi_subset)

    for fasta in odb_generator:
        yield fasta
    for fasta in ncbi_generator:
        yield fasta

def filtered_generator_wrapper(filtered_ids, infile_path, ordered, missing_warning=False):
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
        # filtered = _ordered_filtered_generator(filtered_ids, infile_path)
        filtered = ordered_record_generator(infile_path,filtered_ids)
    else:
        filtered = record_generator(infile_path, filtered_ids)
        # filtered = _filtered_generator(filtered_ids, infile_path)
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
        filtered = filtered_generator_wrapper(filtered_ids,infile_path,ordered)
        SeqIO.write(filtered, outfile_path, "fasta")
    with warnings.catch_warnings(record=True) as w:
        filtered = filtered_generator_wrapper(filtered_ids, infile_path, ordered, True)
    filtered_srs = pd.Series(index=filtered_ids)
    for fasta in filtered:
        filtered_srs[fasta.id] = str(fasta.seq)
    return filtered_srs

###Series, fasta, align_df functions###

def srs_to_fasta(seq_srs, outfile_path):
    # Write records in seq_srs to outfile_path in fasta format
    #Deprecated; use record generator functions to write sequence data (maintains description/ other non-sequence info)
    def record_generator(seq_srs):
        for idx, seq in seq_srs.iteritems():
            record = SeqRecord(Seq(seq, IUPAC.protein), id=idx)
            yield record

    records = record_generator(seq_srs)
    SeqIO.write(records, outfile_path, "fasta")

def fasta_to_srs(fasta_path):
    #Creates series mapping record id to sequence from fasta_path
    with open(fasta_path) as fasta_f:
        fasta_seqs = SeqIO.parse(fasta_f, 'fasta')
        id_seq_map = OrderedDict()
        for fasta in fasta_seqs:
            record_id = fasta.id
            seq = str(fasta.seq)
            id_seq_map[record_id] = seq
        return pd.Series(name="seq", data=id_seq_map)

def align_srs_to_df(align_srs):
    """Returns DataFrame object from series of aligned sequences; columns are 1-indexed positions
    Values are characters in alignment, indexed on record_ids"""
    seq_len = len(align_srs.iloc[0])
    align_df = pd.DataFrame(index=align_srs.index, columns=range(seq_len))
    for idx, seq in align_srs.iteritems():
        align_df.loc[idx, :] = list(seq)
    align_df.columns += 1
    return align_df

def align_fasta_to_df(fasta_path):
    align_srs = fasta_to_srs(fasta_path)
    align_df = align_srs_to_df(align_srs)
    return align_df

def seq_srs_to_align_df(seq_srs, align_in_fpath, align_out_fpath):
    """Transform seq_srs (pandas Series containing sequence texts) to a DataFrame for which each column
    is an alignment position and column. Writes input fasta and output fastas for alignment to align_in_fpath
    and align_out_fpath respectively. Also returns average (non-diagonal) identity distances"""
    srs_to_fasta(seq_srs, align_in_fpath)
    n, ordered_ids, id_dm, align_srs = construct_id_dm(seq_srs, align_in_fpath, align_out_fpath)
    align_df = align_srs_to_df(align_srs)
    return align_df

def align_srs_to_seq_srs(align_srs, outfile_path=None):
    """Return new Series (same index) of sequences with gap characters dropped
    If outfile_path is provided, write un-aligned record seqs to new fasta file"""
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
    tsv_df.drop_duplicates(inplace=True)
    return tsv_df


### Distance Matrix Functions ###
def construct_id_dm(seq_df, seq_fpath, align_outpath="tmp/iddm_align.fasta",
                    ordered=False,aligned=False,kalign_silent=True):
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
    if not aligned:
        # KAlign sequences in filtered_outpath, write to align_outpath
        with open(filtered_fpath,'r') as filtered_f, open(align_outpath,'wt',encoding='utf-8') as align_f:
            args = ['kalign']
            if kalign_silent:
                subprocess.run(args=args, stdin=filtered_f, stdout=align_f, stderr=subprocess.PIPE, text=True)
            else:
                subprocess.run(args=args, stdin=filtered_f, stdout=align_f, text=True)
    else:
        align_outpath = filtered_fpath
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
    #index is a pandas Index object with entries corresponding to the distmat (i.e. lengths and order should be equal)
    #Calculate mean of non-self record distances (diagonal distances generally force-set to 0, so
    #sum functions as intended)
    n = len(distmat)
    avg_dists = np.sum(distmat, axis=1)/(n-1)
    dist_srs = pd.Series(data=avg_dists,index=index,name="dist")
    return dist_srs

def length_srs(fasta_fpath,id_subset=[]):
    """Load sequence and length series corresponding to sequences in fasta_fpath, limited to id_subset if provided

    :param fasta_fpath: File path to fasta of sequences to load length information for
    :param id_subset: if provided, returned series will only contain records in id_subset.
    :return: Series indexed on fasta ids where values are length of record sequences
    """
    fasta_f = open(fasta_fpath)
    fasta_records = SeqIO.parse(fasta_f,'fasta')
    length_dict = {}
    seq_dict = {}
    for fasta in fasta_records:
        fasta_id = fasta.id
        if (len(id_subset) == 0) or \
            (len(id_subset) > 0 and fasta_id in id_subset):
            length_dict[fasta_id] = len(str(fasta.seq))
            seq_dict[fasta_id] = str(fasta.seq)
    lengths = pd.Series(data=length_dict,name='length')
    seqs = pd.Series(data=seq_dict,name='seq')
    fasta_f.close()
    return seqs,lengths
