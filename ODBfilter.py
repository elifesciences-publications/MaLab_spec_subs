from SSerrors import load_errors, SequenceDataError, write_errors
import numpy as np
import pandas as pd
import os
import re
import SSfasta
from IPython.display import display
import warnings

def format_odb_field(field):
    """Remove spaces, commas, and capitalization from alias/ odb fields to search for string matches.
    If field is empty (np.nan), return empty string"""
    if (type(field)) == str:
        field = field.replace(" ", "")
        field = field.replace(",", "")
        field = field.replace("\n", "")
        return field.lower()
    elif type(field) == float and np.isnan(field):
        return ""

def odb_field_to_re(field):
    """Add escape characters for any re special characters present in formatted odb fields."""
    special_chars = ["-","(",")","."]
    field_re = field
    for special_char in special_chars:
        field_re = field_re.replace(special_char,r'\{0}'.format(special_char))
    return field_re


def find_alias_matches(symbol, tsv_df, errors_fpath):
    """Returns a list of the orthodb ids of the reference sequences from an OrthoDB tsv_df and a set containing
    symbol and the GeneCards primary alias for symbol (if it differs from symbol)
    These reference sequences are defined as records with pub_gene_id, og_, or description having a
    text match to either symbol or one of the GeneCards listed aliases for symbol.
    Alias data is fetched from GeneCards automatically and stored in the aliases directory as text file lists.
    :param symbol: Gene symbol from config["IDFilePath]
    :param tsv_df: Unfiltered DataFrame of records from OrthoDB Query tsv for symbol (because function will only be
    called with a valid tsv_df, this function does not do error handling for failed OrthoDB queries)
    :param errors_fpath: File path for error log for run (used to check for failed alias downloads)
    :return am_ids: list of index values from tsv_df for which one of the GeneCards aliases matched the field value
    in tsv_df for pub_gene_id, og_, or description.
    :return exact_matches: list containing accepted gene symbol exact matches for symbol. Contains symbol and optionally
    GeneCards primary alias if different from symbol
    """
    aliases_dir = "aliases_data"
    am_ids = []
    aliases_fpath = "{0}/{1}_aliases.txt".format(aliases_dir,symbol)
    check_errors_file, gc_errors_df = load_errors(errors_fpath,"GeneCardsError")
    if (check_errors_file and symbol in gc_errors_df["gene"].unique()) or not os.path.exists(aliases_fpath):
        #No supplementary GeneCards alias informationl; matches only against symbol
        aliases = [symbol]
        exact_matches = [symbol]
    else:
        # Read previously downloaded alias data from specified dir
        aliases_f = open(aliases_fpath, 'r')
        aliases = [alias.strip() for alias in aliases_f.readlines()]
        gc_name = aliases[0]
        aliases_f.close()
        if gc_name != symbol:
            exact_matches = [symbol,gc_name.upper()]
        else:
            exact_matches = [symbol]
    # Remove spaces, commas, new line chars, and capitalization from alias strings; escape special chars, join into
    #re pattern that can be used with formatted alias fields
    formatted_aliases = [format_odb_field(alias) for alias in aliases]
    alias_REs = [odb_field_to_re(formatted_alias) for formatted_alias in formatted_aliases]
    aliases_pat = "({0})".format("|".join(alias_REs))
    # Search fields in search_fields for matches to the alias strings provided by GeneCards
    # Iterate tsv_df rows, save all reference ids which have matches
    search_fields = ["pub_gene_id", "og_name", "description"]

    for idx, row in tsv_df.iterrows():
        #         for field in search_fields:
        # Current behavior: exact string matches in formatted pub_gene_id, og_name, or description only to one of the
        #GeneCards aliases. Second function for generating exact matches to exact_matches (IDFile given gene symbol or
        #GeneCards primary alias only).
        # TODO: Add in partial alias string matching. Difficulties with distinguishing symbols

        for search_field in search_fields:
            field_val = row[search_field]
            if type(field_val)==float and np.isnan(field_val):
                #Ignore empty values in tsv_df rows
                continue
            formatted_field = format_odb_field(field_val)
            if idx not in am_ids and re.search(aliases_pat,formatted_field):
                am_ids.append(idx)
    # exact_matches contains only IDFile stored gene symbols or GeneCards primary aliases
    if len(am_ids) == 0:
        msg = 'No alias matched sequences could be found in OrthoDB input for gene symbol {0}'.format(symbol)
        raise SequenceDataError(0,msg)
    return am_ids, exact_matches

def exact_match_df(unfiltered_df,exact_matches):
    """From an unfiltered_df of reference sequences, returns a DataFrame of all entries which have pub_gene_id
    matching symbols provided in exact_matches.

    :param unfiltered_df:
    :param exact_matches:
    :return final_df: DataFrame containing all records per species in organism_taxid corresponding to
    an exact match of pub_gene_id to exact_matches. 
    """
    #Upper matches: exact matches reg. exps matching entries in exact_matches followed only by end of string or non
    #alpha-numeric ([^\w]) characters to prevent partial matching against different gene symbols
    upper_matches = [match.upper() for match in exact_matches]
    # upper_matches = [match+"$|"+match+"[^\w]" for match in upper_matches]
    upper_matches = [r'{0}$|{0}[^\w]'.format(match) for match in upper_matches]
    pat = "|".join(upper_matches)
    pg_id_df = unfiltered_df.loc[unfiltered_df["pub_gene_id"].str.upper().str.contains(pat)]
    #Filtering for any species with more than one sequence
    final_df = pd.DataFrame(columns=pg_id_df.columns)
    for unique_tax in pg_id_df["organism_taxid"].unique():
        spec_df = pg_id_df.loc[pg_id_df["organism_taxid"] == unique_tax]
        # if len(spec_df) > 1:
        #     min_dist_idx = spec_df["dist"].values.argmin()
        #     min_dist_row = spec_df.iloc[min_dist_idx,:]
        #     final_df = final_df.append(min_dist_row)
        # else:
        #     final_df = final_df.append(spec_df)
        final_df = final_df.append(spec_df)
    return final_df

def __parse_manual_selection_input(gene_symbol,selection_df,display_df,manual_selections_fpath):
    print("Matched records, choose representative input sequence from below.")
    display(display_df)

    #Initialize cached manual selections table file, check for previous selections
    if os.path.exists(manual_selections_fpath):
        ms_df = pd.read_csv(manual_selections_fpath,sep='\t',index_col='gene_symbol')
        if gene_symbol in ms_df.index:
            record_id = ms_df.loc[gene_symbol,'record_id']
            selection_row = selection_df.loc[record_id,:]
            print("Using cached record selection: record_id {0}".format(record_id))
            print("To clear selections, either delete corresponding row in file at {0} ".format(manual_selections_fpath))
            return selection_row
    else:
        ms_df = pd.DataFrame(columns=['record_id'],dtype=str)
        ms_df.index.name = 'gene_symbol'

    while True:
        try:
            input_idx = input("Enter 0-indexed position of representative sequence for "
                              "analysis (ie 0 for first row, 1 for second)")
            int_idx = int(input_idx)
            if int_idx < 0:
                #Special case because iloc will accept negative integer indices
                raise ValueError
            selection_row = selection_df.iloc[int_idx, :]
        except (IndexError, ValueError) as e:
            print("Couldn't parse input, enter a number between 0 and {0}".format(len(selection_df)-1))
            continue
        else:
            break
    record_id = selection_row.name
    ms_df.loc[gene_symbol,'record_id'] = record_id
    ms_df.to_csv(manual_selections_fpath,sep='\t')
    return selection_row

def min_dist_spec_record(distmat,dm_record_ids,spec_record_ids,against_record_ids, record_df):
    """Calculates the average distance of every record containing spec_taxid against accepted records, then
    returns the row from ref_df corresponding to the record with lowest average distance.

    :param distmat: np.ndarray n x n distance matrix
    :param dm_record_ids: iterable of record_ids corresponding to rows of distmat
    :param spec_record_ids: list of records from specific species, subset of dm_record_ids
    :param against_record_ids: record ids against which record distances will be calculated, subset of dm_record_ids
    :param record_df: Sequence dataframe with correspinding records to distmat
    :return: md_row: row from ref_df with minimum average distance to against_record_id records,
    :return min_dist: minimum average distance value
    """
    #Issue warnings for record IDs from spec_record_ids and against_record_ids
    for rid in spec_record_ids:
        if rid not in dm_record_ids:
            msg = "Species record ID {0} not in dm_record_ids.".format(rid)
            warnings.warn(msg)
    for rid in against_record_ids:
        if rid not in dm_record_ids:
            msg = "Against record ID {0} not in dm_record_ids.".format(rid)
            warnings.warn(msg)
    #Identify index positions in distmat for record IDs for species/ against records
    spec_records = [(i,id_) for i,id_ in enumerate(dm_record_ids) if id_ in spec_record_ids]
    spec_dm_idxs = [t[0] for t in spec_records]
    against_records = [(i,id_) for i,id_ in enumerate(dm_record_ids) if id_ in against_record_ids]
    accepted_dm_idxs = [t[0] for t in against_records]
    #sub_dm: distmat submatrix of spec_record columns and against_record rows
    spec_dm = distmat[:,spec_dm_idxs]
    sub_dm = spec_dm[accepted_dm_idxs,:]
    #Flatten into [n_spec_records] array, identify min_idx and return corresponding spec_record row and distance
    if len(sub_dm) > 1:
        avg_dist = sub_dm.mean(axis=0)
    else:
        avg_dist = sub_dm[0]
    min_idx = np.argmin(avg_dist)
    min_dist_id = spec_records[min_idx][1]
    min_dist = avg_dist[min_idx]
    md_row = record_df.loc[min_dist_id,:]
    return md_row, min_dist

def select_known_species_records(gene_symbol,em_df, am_df, ks_taxids, ks_refseqs_fpath,
                                 manual_selections_fpath = 'tmp/manual_record_selections.tsv'):
    """Return a dataframe of at most one record per species in ks_taxids of representative sequences for species in
    ks_taxids. ks_taxids will by default be set to include well-annotated species (human/mouse) and the test species
    from the config file (by default 13LGS).
    Representative sequences are interpreted as being the only avialable record from a species or the sequence
    most similar (max identity) to other ks_taxids rep. sequences if a species has multiple records. Sequences which are
    not present in am_df (ie no alias matches) will not be returned even if present for a species in the unfiltered input.
    Sequences will be searched for first in em_df but in am_df if no records for that species are in em_df. If em_df or
    am_df have exactly one record per species, those records will be assumed to be correct and no distance calculations
    or filtering will be done. Returned DataFrame can possibly have no records present for test_species, in which case
    the function calling select_known_species_records must check for test_species record presence.

    :param gene_symbol: String HGNC identifier for gene of input data
    :param (DataFrame) em_df: DataFrame of exact symbol matched records
    :param (DataFrame) am_df: DataFrame of alias matched records
    :param (str) ts_taxid: Taxonomy ID for test species
    :param (array-like) ks_taxids: Taxonomy IDs for well-annotated species and test species (human/ mouse/ 13LGS)
    :param ks_refseqs_fpath: Fasta file path for alias-match filtered OrthoDB input
    :return:
    """
    #Filter em_df and am_df down to taxonomy IDs in ks_taxids
    ksr_am_df, ksr_em_df = am_df.loc[am_df["organism_taxid"].isin(ks_taxids),:],\
                           em_df.loc[em_df["organism_taxid"].isin(ks_taxids), :]
    am_taxid_uniques,em_taxid_uniques = ksr_am_df["organism_taxid"].unique(),ksr_em_df["organism_taxid"].unique()
    final_ksr_df = pd.DataFrame(columns=ksr_em_df.columns)
    #If exactly 3 records from exactly 3 species, reorder in order of ks_taxids and return. Saves extra alignment steps
    if len(ksr_em_df) == len(ks_taxids) and len(em_taxid_uniques) == len(ks_taxids):
        for tax_id in ks_taxids:
            row = ksr_em_df.loc[ksr_em_df["organism_taxid"] == tax_id, :]
            final_ksr_df = final_ksr_df.append(row)
        return final_ksr_df
    elif len(ksr_am_df) == len(ks_taxids) and len(am_taxid_uniques) == len(ks_taxids):
        for tax_id in ks_taxids:
            row = ksr_am_df.loc[ksr_am_df["organism_taxid"] == tax_id, :]
            final_ksr_df = final_ksr_df.append(row)
        return final_ksr_df
    #Set selection df to be the smallest available set for which at least one record present from ks_taxids species
    if ksr_em_df.empty:
        if ksr_am_df.empty:
            raise SequenceDataError(1, "No GeneCards alias matched sequence records for human/mouse/test species")
        else:
            selection_df = ksr_am_df
    else:
        selection_df = ksr_em_df
    #Populate single_avail_ksr with records if there is only one record in selection_df from that tax_id
    single_avail_ksr = pd.DataFrame(columns=ksr_em_df.columns)
    for tax_id in ks_taxids:
        tax_em_df = selection_df.loc[selection_df['organism_taxid']==tax_id,:]
        if len(tax_em_df) == 1:
            single_avail_ksr = single_avail_ksr.append(tax_em_df)
    if single_avail_ksr.empty:
        #If no species have single record, take manual input (or read from cached selections if previously entered),
        #use selected record as seed input for determining best records from other species.
        selection_fapath = 'tmp/filtered_selection_intput.fasta'
        SSfasta.filter_fasta_infile(selection_df.index,ks_refseqs_fpath,selection_fapath)
        display_df = selection_df.copy().drop(columns=['pub_og_id','og_name','level_taxid'])
        display_df.loc[:,'seq'] = SSfasta.fasta_to_srs(selection_fapath)
        selection_row = __parse_manual_selection_input(gene_symbol,selection_df,display_df,manual_selections_fpath)
        single_avail_ksr = single_avail_ksr.append(selection_row)
    sa_record_ids = single_avail_ksr.index
    sa_taxid_uniques = single_avail_ksr['organism_taxid'].unique()

    am_id_dm, am_align_srs = SSfasta.construct_id_dm(ksr_am_df, ks_refseqs_fpath,ordered=False)

    for ks_id in ks_taxids:
        if ks_id not in sa_taxid_uniques:
            #Use em_df or am_df depending on if ks_id is present
            if ks_id in em_taxid_uniques:
                spec_record_ids = ksr_em_df.loc[ksr_em_df['organism_taxid'] == ks_id, :].index
            elif ks_id in ks_id in am_taxid_uniques:
                spec_record_ids = ksr_am_df.loc[ksr_am_df['organism_taxid'] == ks_id, :].index
            else:
                #If no records for taxid in either em or am dfs, skip ksr selection
                continue
            # Maximum identity = minimum id_dm value based on AlignIO implementation
            md_row, min_dist = min_dist_spec_record(am_id_dm, am_align_srs.index, spec_record_ids,
                                                    sa_record_ids, ksr_am_df)
            final_ksr_df = final_ksr_df.append(md_row)
        else:
            sa_row = single_avail_ksr.loc[single_avail_ksr['organism_taxid'] == ks_id, :]
            final_ksr_df = final_ksr_df.append(sa_row)
    return final_ksr_df


def select_outgrup_records(em_df, am_df, ks_taxids,final_ksr_df, seqs_fpath):
    """Select records for remaining OrthoDB outgroup species in analysis that are not in ks_taxids.

    Selection is based on maximum identity to accepted records in final_ksr_df (ie accepted human/mouse/13LGS); best
    (max identity) alias-matched records for a species will still be dropped if they do not meet an identity threshold,
    currently set to 1.5 * (average identity of sequences in final_lsr_df against each other).

    :param (DataFrame) em_df: DataFrame of records with exact pub_gene_id match to accepted gene symbols
    :param (DataFrame) am_df: DataFrame of records matching aliases for symbol
    :param (array-like) ks_taxids: list of taxonomy ids to use as accepted species records
    :param final_ksr_df: DataFrame of accepted records from well-annotated species (max one per species).
    :param seqs_fpath: Fasta file path containing at least all records in am_df. Can be safely set to unfiltered
    fasta input, records will be automatically filtered down appropriately using am_df.
    :return final_df: DataFrame containing all selected records (final_ksr_df records first)
    :return dist_srs: Series indexed on record_id and containing average distance information against rest of input set.
    """
    am_non_ksr_taxids = [tax_id for tax_id in am_df["organism_taxid"].unique() if tax_id not in ks_taxids]

    # Distance calculations for final set of known species records - check internal identity values
    # Set identity threshold - other species sequences above this value will not be included
    am_dm_fpath = "tmp/am_dm_ka.fasta"
    am_id_dm,am_align_srs = SSfasta.construct_id_dm(am_df,seqs_fpath,am_dm_fpath)
    am_record_idx = am_align_srs.index
    ksr_record_idx = final_ksr_df.index
    ksr_pos = [am_record_idx.get_loc(record_id) for record_id in ksr_record_idx]
    n_ksr = len(ksr_pos)
    #ksr_sub_dm: n_ksr x n_ksr distance matrix consisting of values of ksr records against each other
    ksr_sub_dm = am_id_dm[:,ksr_pos]
    ksr_sub_dm = ksr_sub_dm[ksr_pos,:]
    # Ignore diagonal (0 values for record identity against itself); set identity threshold which is used to exclude
    #dissimilar records
    non_diagonal_avg = ksr_sub_dm.sum(axis=0) / (n_ksr - 1)
    identity_threshold = np.mean(non_diagonal_avg) * 1.5
    final_df = final_ksr_df.copy()
    #Add in min dist records for other outgroup species (am_non_ksr_taxids) to final_df if they meet identity threshold
    for taxid in am_non_ksr_taxids:
        tax_records = am_record_idx[am_record_idx.str.contains(taxid)]
        tax_pos = [am_record_idx.get_loc(record_id) for record_id in tax_records]
        md_row, md = min_dist_spec_record(am_id_dm,am_record_idx,tax_records,ksr_record_idx,am_df)
        if md <= identity_threshold:
            final_df = final_df.append(md_row)
        else:
            print("Min dist record for tax_id {0} does not meet distance threshold {1}".format(taxid,identity_threshold))

    #Reorder and filter distmat down to final_df records order, calculate non-diagonal avg distances
    dm_pos = [am_record_idx.get_loc(final_idx) for final_idx in final_df.index]
    final_ordered_dm = am_id_dm[dm_pos,:]
    final_ordered_dm = final_ordered_dm[:,dm_pos]
    dist_srs = SSfasta.avg_dist_srs(final_df.index,final_ordered_dm)

    return final_df, dist_srs

def final_ksr_df_QC(gene_symbol,matches,final_ksr_df,ks_taxids,ts_taxid,
                    seq_qc_fpath,seq_fpath,length_warning=False):
    """Log quality checks on final accepted record sequences for species in ks_taxids.

    Entries which fail QC checks will still be run in analysis, but the quality checks for sequence consistency logged
    may serve as warnings by which the user can elect to remove results for those gene symbols. The only QC which will
    prevent analysis is absence of test species tax_id in the final accepted records dataframe.

    :param gene_symbol: Gene symbol for which data is being used
    :param matches: List of acceptable exact symbol matches for gene_symbol
    :param seq_qc_fpath: Log file for quality checks
    :param seq_fpath: Path to fasta with sequences
    :param final_ksr_df: DataFrame returned by select_known_species_records
    :param ks_taxids: list of tax_ids expected to be in final_ksr_df
    :return: N/A
    """
    if len(final_ksr_df) < len(ks_taxids):
        if ts_taxid not in final_ksr_df["organism_taxid"].unique():
            msg = "No alias matched sequence could be found for test species (taxid: {0})".format(ts_taxid)
            raise SequenceDataError(2,msg)
        for tax_id in ks_taxids:
            if tax_id not in final_ksr_df["organism_taxid"].unique():
                message_txt = "No reference sequence for tax_id: {0}".format(tax_id)
                write_ref_seq_QC(seq_qc_fpath,gene_symbol,message_txt)
    if length_warning:
        # length_srs = final_ksr_df["length"]
        length_srs = SSfasta.length_srs(seq_fpath,final_ksr_df.index)
        median_len = length_srs.median()
        for record_id in final_ksr_df.index:
            id_len = length_srs[record_id]
            if (np.abs(id_len-median_len)/median_len) >= 0.1:
                message_txt = "Record_id {0} has length {1} which is greater than 10% different from the median ({2})".format(record_id,id_len,median_len)
                write_ref_seq_QC(seq_qc_fpath,gene_symbol,message_txt)
    upper_matches = [match.upper() for match in matches]
    upper_matches = [match+"$|"+match+"[;]" for match in upper_matches]
    pat = "|".join(upper_matches)
    for record_id,pgid in final_ksr_df["pub_gene_id"].iteritems():
        if not re.search(pat,pgid.upper()):
            message_txt = "Record_id {0} has pub_gene_id {1} which doesn't match gene_symbol ({2})".format(record_id,pgid,gene_symbol)
            write_ref_seq_QC(seq_qc_fpath,gene_symbol,message_txt)

def write_ref_seq_QC(seq_qc_fpath,gene_symbol,message):
    """Writes warning messages to seq_qc_fpath and prints (will not write duplicate message entries)

    :param seq_qc_fpath: Log file for quality check messages
    :param gene_symbol:
    :param message:
    :return:
    """
    if not os.path.exists(seq_qc_fpath):
        seq_qc_df = pd.DataFrame(columns=['gene_symbol','message'])
    else:
        seq_qc_df = pd.read_csv(seq_qc_fpath,sep='\t',index_col=0)
        if gene_symbol in seq_qc_df['gene_symbol'].unique():
            symbol_df = seq_qc_df[seq_qc_df['gene_symbol']==gene_symbol]
            for idx,row in symbol_df.iterrows():
                row_symb, row_msg = row.values
                if row_msg == message:
                    print("{0}\t{1}".format(row_symb,row_msg))
                    return
    seq_qc_df = seq_qc_df.append(pd.Series({'gene_symbol':gene_symbol,'message':message}),ignore_index=True)
    seq_qc_df.to_csv(seq_qc_fpath,sep='\t')
    print("{0}\t{1}".format(gene_symbol, message))


def process_input(symbol,config,tax_subset):
    """Return final

    :param symbol:
    :param config:
    :param tax_subset:
    :return:
    """

    run_name, test_tid = config['RunName'], config['ODBTestTaxID']
    raw_tsv_fpath,raw_fa_fpath = "{0}/input/ODB/{1}.tsv".format(run_name,symbol),\
                                 "{0}/input/ODB/{1}.fasta".format(run_name,symbol)
    seq_qc_fpath = '{0}/summary/accepted_record_QC.tsv'.format(run_name)
    errors_fpath = '{0}/summary/error.tsv'.format(run_name)
    ks_taxids = ['10090_0', '43179_0', '9606_0']
    # tsv_inpath = "cDNAscreen_041020/input/ODB/{0}.tsv".format(symbol)
    unfiltered_tsv = SSfasta.load_tsv_table(raw_tsv_fpath, tax_subset=tax_subset)
    #Filter by alias matches, exact pub_gene_id matches
    try:
        am_ids, exact_matches = find_alias_matches(symbol, unfiltered_tsv, errors_fpath)
        am_df = unfiltered_tsv.loc[am_ids]
        em_df = exact_match_df(unfiltered_tsv, exact_matches)
        final_ksr_df = select_known_species_records(symbol, em_df, am_df, ks_taxids, raw_fa_fpath)
        final_ksr_df_QC(symbol,exact_matches,final_ksr_df,ks_taxids,test_tid,seq_qc_fpath,raw_fa_fpath)
        final_input_df, dist_srs = select_outgrup_records(em_df, am_df, ks_taxids, final_ksr_df, raw_fa_fpath)
    except SequenceDataError as sde:
        #Log errors, raise error for handling in calling function
        write_errors(errors_fpath,symbol,sde)
        raise sde
    return final_input_df


def main():
    pass

if __name__ == "__main__":
    main()