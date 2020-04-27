from SSerrors import load_errors, SequenceDataError
import numpy as np
import pandas as pd
import os
import re
import SSfasta
from IPython.display import display
import warnings
from Bio import SeqIO,Seq
import SSdirectory

from ODBfilter import min_dist_spec_record, select_known_species_records

def ordered_record_generator(fpath,ordered_ids):
    """Generates Seq records from fpath, limited to ids if provided. Generator order is order in ordered_ids.
    :param fpath: Fasta file path
    :param ordered_ids: array-like containing ordered record ids
    :return:
    """
    for id in ordered_ids:
        with open(fpath) as fasta_f:
            fastas = SeqIO.parse(fasta_f,'fasta')
            for fasta in fastas:
                if fasta.id == id:
                    yield fasta
                    break

def record_generator(fpath,ids=[]):
    """Generates Seq records from fpath, limited to ids if provided. Ordered as in fpath
    :param fpath: Fasta file path
    :param ids: If provided, only yield records corresponding to record ids provided. Default value causes all records
    to be present in generator
    :return: generator object
    """
    with open(fpath) as fasta_f:
        fastas = SeqIO.parse(fasta_f, 'fasta')
        for fasta in fastas:
            if (len(ids)>0 and fasta.id in ids) or \
                    (len(ids)==0):
                yield fasta

def ODB_NCBI_generator(ODB_fpath,NCBI_fpath,odb_subset=[],ncbi_subset=[],ordered=False):
    """

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

def load_NCBI_fasta_df(NCBI_fasta_fpath,taxid_dict):
    """Reads NCBI fasta into DataFrame, extracting available fields into appropritate columns

    :param NCBI_fasta_fpath: File path to NCBI fasta
    :param taxid_dict: maps species names to NCBI_taxids
    :return: DataFrame populated with record information from NCBI fasta 
    """
    ncbi_fastas = SeqIO.parse(NCBI_fasta_fpath, 'fasta')
    ncbi_df = pd.DataFrame(columns=["organism_taxid", "organism_name", "description", "length", "seq"])
    for fasta in ncbi_fastas:
        row_dict = {}
        desc_remaining = re.search("{0}(.*)".format(fasta.id), fasta.description).groups()[0]
        if desc_remaining:
            #If standard format description string, will extract species name and description.
            desc_remaining = desc_remaining.strip()
            organism_name = re.search(r"\[(\w+\s\w+(\s\w+)?)\]$", desc_remaining).groups()[0].strip()
            organism_taxid = taxid_dict[organism_name]
            row_dict['organism_name'],row_dict['organism_taxid'] = organism_name,organism_taxid
            row_dict['description']  = desc_remaining
        row_dict['seq'],row_dict['length']= str(fasta.seq),len(str(fasta.seq))
        f_row = pd.Series(data=row_dict,name=fasta.id)
        ncbi_df = ncbi_df.append(f_row)
    return ncbi_df

def select_NCBI_record(ODB_fasta_fpath,NCBI_fasta_fpath,taxid_dict,ODB_final_input_df,compare_taxids):
    """Selects best NCBI record from NCBI fasta fpath by max identity to the OrthoDB records represented by compare_taxids.

    :param ODB_fasta_fpath: Fasta path for ODB records
    :param NCBI_fasta_fpath: Fasta path to NCBI records (should only contain records from one species)
    :param taxid_dict: Maps species names to taxids, used by load_NCBI_fasta_df
    :param ODB_final_input_df: DataFrame containing accepted records from OrthoDB data, returned from
    ODBfilter.process_input
    :param (collection) compare_taxids: tax_ids against which distance should be calculated to determine minimum
    distance NCBI record
    :return: combined_df, DataFrame containing rows from ODB_final_input and the minimu, distance row from
    NCBI_fasta_fpath
    """
    
    ncbi_df = load_NCBI_fasta_df(NCBI_fasta_fpath,taxid_dict)
    if len(ncbi_df) > 1:
        #Align all unfiltered NCBI records against ODB_final_input records
        combined_unaln_fpath,combined_aln_fpath = "tmp/ODB_NCBI_unaln.fasta","tmp/ODB_NCBI_aln.fasta"
        unaln_generator = ODB_NCBI_generator(ODB_fasta_fpath,NCBI_fasta_fpath,odb_subset=ODB_final_input_df.index)
        SeqIO.write(unaln_generator, combined_unaln_fpath, "fasta")
        combined_df = ODB_final_input_df.append(ncbi_df,sort=False)
        # display(combined_df)
        id_dm, align_srs = SSfasta.construct_id_dm(combined_df,combined_unaln_fpath,align_outpath=combined_aln_fpath)
        spec_record_ids= ncbi_df.index
        compare_record_ids = ODB_final_input_df.loc[ODB_final_input_df['organism_taxid'].isin(compare_taxids)].index
        md_row,min_dist = min_dist_spec_record(id_dm,align_srs.index,spec_record_ids,compare_record_ids,combined_df)
        final_combined_df = ODB_final_input_df.append(md_row,sort=False)
        return final_combined_df
    else:
        combined_df = ODB_final_input_df.append(ncbi_df,sort=False)
        return combined_df

def annotate_source_and_filter(config,symbol,combined_df,am_df,em_df,
                               manual_selections_fpath="tmp/manual_record_selections.tsv",
                               drop_cols=['pub_og_id','level_taxid']):
    """Populates source_db and filter_type information into combined_df, reorders cols and returns modified DataFrame.

    :param combined_df: DataFrame containing final record set from OrthoDB and NCBI
    :param am_df: alias matched OrthoDB records
    :param em_df: exact symbol matched OrthoDB records
    :param manual_selections_fpath: tsv file containing records which were manually selected. If file exists at provided
    path, changes selection_type value in corresponding record row
    :param (array-like) drop_cols: If provided, drops columns from processed DataFrame. Used currently to filter
    some redundant OrthoDB tsv information.
    :return: Returns edited combined_df
    """
    #Leave original dataframe untouched
    combined_df = combined_df.copy()
    run_name, ncbi_taxid = config['RUN']['RunName'], config['NCBI']['NCBITaxID']
    selection_col = 'selection_type'

    odb_records_idx = combined_df.index[combined_df.index.isin(am_df.index)]
    ncbi_records_idx = combined_df.index[~combined_df.index.isin(am_df.index)]
    combined_df.loc[odb_records_idx, 'db_source'] = "OrthoDB"
    combined_df.loc[ncbi_records_idx, 'db_source'] = "NCBI"

    unfiltered_NCBI_fasta = "{0}/input/NCBI/{1}/{2}.fasta".format(run_name, ncbi_taxid, symbol)
    unf_ncbi_srs = SSfasta.fasta_to_srs(unfiltered_NCBI_fasta)
    if len(unf_ncbi_srs) > 1:
        ncbi_filt = "NCBI min dist"
    else:
        ncbi_filt = "NCBI single record"
    combined_df.loc[ncbi_records_idx, selection_col] = ncbi_filt

    for record_id,row in combined_df.loc[odb_records_idx,:].iterrows():
        taxid = row['organism_taxid']
        em_taxid_df = em_df.loc[em_df['organism_taxid']==taxid,:]
        am_taxid_df = am_df.loc[am_df['organism_taxid']==taxid,:]
        if len(em_taxid_df) == 0:
            if len(am_taxid_df) == 1:
                combined_df.loc[record_id,selection_col] = "alias match single record"
            else:
                combined_df.loc[record_id, selection_col] = "alias match min dist"
        elif len(em_taxid_df) == 1:
            combined_df.loc[record_id, selection_col] = "symbol match single record"
        else:
            combined_df.loc[record_id,selection_col] = "symbol match min dist"
    #Read manual selection information, fix selection_type if exists for symbol
    if os.path.exists(manual_selections_fpath):
        manual_selections_df = pd.read_csv(manual_selections_fpath,sep='\t',dtype=str,index_col='gene_symbol')
        if symbol in manual_selections_df.index:
            record_id = manual_selections_df.loc[symbol,'record_id']
            combined_df.loc[record_id,selection_col] = 'manual selection'
    #Drop columns from drop_cols
    if drop_cols:
        combined_df.drop(columns=drop_cols,inplace=True)
    #Reorder columns so seq is last.
    reordered_cols = combined_df.columns.tolist()
    reorder_labels = ['db_source','selection_type','length','dist','seq']
    for label in reorder_labels:
        label_pos = reordered_cols.index(label)
        reordered_cols.pop(label_pos)
        reordered_cols.append(label)
    combined_df = combined_df[reordered_cols]
    return combined_df



def combined_records_processing(config,am_df,em_df,combined_df,
                                symbol, odb_fasta="", ncbi_fasta="",out_unaln_fasta="",out_aln_fasta="",
                                out_tsv_fpath=""):
    """Adds in internal record distance information and source DataBase annotations for final dataset. Writes final
    input dataset (both OrthoDB and NCBI records) to 1) unaligned fasta 2) aligned fasta and 3) a records table
    corresponding to the record modified record DataFrame

    :param config: configparser object
    :param am_df alias_match DataFrame (see ODBfilter)
    :param em_df: exact match DataFrame
    :param combined_df: ODB and NCBI combined record DataFrame as returned by select_NCBI_record
    :param symbol: Gene symbol
    :param odb_fasta,ncbi_fasta: If provided, will use as sources for records to write final dataset sequences from. If
    not provided, defaults to unfiltered OrthoDB and NCBI fasta paths given run_name and symbol
    :param out_unaln_fasta,out_aln_fasta: If provided, will write unaligned/aligned final record Sequence set to
    these paths. If not provided, uses default output directory path for symbol.
    :return: processed_df: DataFrame modified to include record distance, db_source, and filter_type
    """
    run_name,ncbi_taxid = config['RUN']['RunName'],config['NCBI']['NCBITaxID']
    if not odb_fasta:
        odb_fasta = "{0}/input/ODB/{1}.fasta".format(run_name,symbol)
    if not ncbi_fasta:
        ncbi_fasta = "{0}/input/NCBI/{1}/{2}.fasta".format(run_name,ncbi_taxid,symbol)
    if not out_unaln_fasta:
        out_unaln_fasta = "{0}/output/{1}/{1}.fasta".format(run_name,symbol)
    if not out_aln_fasta:
        out_aln_fasta = "{0}/output/{1}/{1}_msa.fasta".format(run_name, symbol)
    if not out_tsv_fpath:
        out_tsv_fpath = "{0}/output/{1}/{1}_records.tsv".format(run_name, symbol)

    SSdirectory.create_directory("{0}/output/{1}".format(run_name,symbol))
    manual_selections_fpath = "{0}/manual_record_selections.tsv".format(run_name)

    #Final unaligned and aligned Fasta writing
    odb_records_idx = combined_df.index[combined_df.index.isin(am_df.index)]
    ncbi_records_idx = combined_df.index[~combined_df.index.isin(am_df.index)]
    combined_records = ODB_NCBI_generator(odb_fasta,ncbi_fasta,odb_subset=odb_records_idx,
                                          ncbi_subset=ncbi_records_idx,ordered=True)
    SeqIO.write(combined_records, out_unaln_fasta, 'fasta')
    #Internal distance calculation
    id_dm, aln_srs = SSfasta.construct_id_dm(combined_df,out_unaln_fasta,align_outpath=out_aln_fasta)
    dist_srs = SSfasta.avg_dist_srs(combined_df.index,id_dm)
    combined_processed = combined_df.copy()
    combined_processed.loc[:,'dist'] = dist_srs
    #Add source_db and selection_type information, drop redundant columns (default is level_taxid, pub_og_id).
    combined_processed = annotate_source_and_filter(config,symbol,combined_processed,am_df,em_df,manual_selections_fpath)
    #write records table to file
    combined_processed.to_csv(out_tsv_fpath,sep='\t')
    return combined_processed
