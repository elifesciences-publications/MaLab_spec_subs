from SSerrors import load_errors, SequenceDataError
import numpy as np
import pandas as pd
import os
import re
import SSfasta
from IPython.display import display
import warnings
from Bio import SeqIO,Seq

from ODBfilter import min_dist_spec_record, select_known_species_records

def ordered_record_generator(fpath,ordered_ids):
    """Generates Seq records from fpath, limited to ids if provided. Generator order is order in ordered_ids.
    :param fpath: Fasta file path
    :param ordered_ids: array-like containing ordered record ids
    :return:
    """
    with open(fpath) as fasta_f:
        fastas = SeqIO.parse(fasta_f,'fasta')
        for id in ordered_ids:
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
    """

    :param ODB_fasta_fpath: Fasta path for ODB records
    :param NCBI_fasta_fpath: Fasta path to NCBI records (should only contain records from one species)
    :param taxid_dict: Maps species names to taxids, used by load_NCBI_fasta_df
    :param ODB_final_input_df:
    :param compare_taxids:
    :return:
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
