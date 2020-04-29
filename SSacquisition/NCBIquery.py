#NCBIquery.py - Handles ortholog Gene ID mapping (selenium) and Protein Record Downloads (entrez API)
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

import os
import pandas as pd
from SSutility.SSerrors import RecordDataError,NCBIQueryError, write_errors, print_errors, load_errors
import re
from selenium import webdriver
from selenium.common.exceptions import TimeoutException,NoSuchElementException
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException,NoSuchElementException
from xml.etree.ElementTree import ElementTree as ET
import numpy as np

def headless_driver():
    #Headless driver with default options
    from selenium.webdriver.chrome.options import Options
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    WINDOW_SIZE = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % WINDOW_SIZE)
    driver = webdriver.Chrome(options=chrome_options)
    return driver


def download_AGS_data(gene_id_df, config):
    """Uses NCBI Rest API to acquire sequence data for NCBI gene IDs listed in gene_id_df[human_gene_id]

    :param gene_id_df: DataFrame object with required columns gene_symbol and human_gene_id
    :param config: configparser object created from config/config.txt
    :param errors_tsv_fpath: filepath to errors table (see SSerrors.write_errors for details)
    :return: ags_geneID_df: DataFrame containing all columns originally in gene_id_df as well as columns
    specified in the config file by NCBIGeneIDField and NCBIProteinIDField with values populated by map_AGS_geneIDs
    and download_NCBI_records respectively.
    """
    from SSutility.SSdirectory import create_directory
    run_config, ncbi_config = config['RUN'],config['NCBI']
    run_name, errors_fname = config['RUN']['RunName'], config['RUN']['ErrorsFileName']
    errors_fpath = "{0}/{1}".format(run_name, errors_fname)
    NCBI_taxid,NCBI_spec_name = ncbi_config["NCBITaxID"], ncbi_config['NCBITaxName']
    gene_field_name, protein_field_name = ncbi_config['NCBIGeneIDField'],ncbi_config['NCBIProteinIDField']
    NCBI_API_key = ncbi_config.get("NCBIAPIKey","")
    NCBI_input_dir = "{0}/input/NCBI/{1}".format(run_name,NCBI_taxid)
    create_directory(NCBI_input_dir)

    filled_outpath = "{0}/summary/cDNAscreen_geneIDs_complete.tsv".format(run_name)

    tax_dict = {'gid_column':gene_field_name,'pid_column':protein_field_name,
                'spec_name':NCBI_spec_name,'taxid':NCBI_taxid}
    mapped_id_df = map_AGS_geneIDs(gene_id_df, filled_outpath, errors_fpath,tax_dict)
    ags_mapped_id_df = download_NCBI_records(mapped_id_df, NCBI_input_dir,tax_dict,NCBI_API_key,
                                             pid_outpath=filled_outpath)
    return ags_mapped_id_df

def single_NCBI_gid_query(driver,id_df,idx,hgid,tax_dict,
                       id_out_fpath="",initial_timeout=1):
    """For a human NCBI gene id, searches for ortholog gene ids corresponding to all taxids in taxid_columns.
    NOTE: Edited from spec_subs version of function to support multiple taxids per gene_id query

    :param driver (selenium.webdriver): Selenium webdriver object, initialized pre-iteration with headless option and implicit wait
    :param (DataFrame) id_df: DataFrame object where ortholog gene IDs will be stored
    :param (idx) idx: for row in id_df in which to store the ortholog gene ID
    :param (str) hgid: NCBI human Gene ID for given gene symbol
    :param (dict) tax_dict: Dictionary mapping 'gid_column' to a column label for id_df and 'spec_name' to scientific
    species name for which NCBI ortholog GID is being queried
    :param (str) id_out_fpath: if provided, will write updated id_df as tsv to provided path.
    :param (int) initial_timeout: Initial timeout length for loading NCBI Gene orthologs page. Default 1

    :return: None. Edits id_df directly, or raises NCBIQueryError to log entries with no ortholog data for any species
    """

    #Mammalian taxonomy level: 40674 - ensure this is a valid taxonomy identifier and encompasses all NCBI species
    # in analysis or all requests will fail to produce ortholog data
    taxonomy_level = '40674'
    req_url = "https://www.ncbi.nlm.nih.gov/gene/{0}/ortholog/?scope={1}".format(hgid,taxonomy_level)
    driver.get(req_url)
    driver_timeout = initial_timeout
    out_col = tax_dict['gid_column']
    spec_name = tax_dict['spec_name']
    gene_xpath = "//tr/td/label/*[contains(text(), '{0}')]".format(spec_name) + \
                 "/parent::label/parent::td/parent::tr/td[@class=' fld-gene']/a"
    try:
        spec_elem = WebDriverWait(driver,driver_timeout).until(
            EC.presence_of_element_located((By.XPATH,gene_xpath))
        )
    except TimeoutException:
        id_df.loc[idx, out_col] = np.nan
        #Check for generic ortholog table row elements (ie presence of any orthologs at all)
        generic_tr_xpath = "//tr/td/label"
        try:
            driver.find_element_by_xpath(generic_tr_xpath)
        except NoSuchElementException:
            # No orthologs at all in table, raise NCBIQueryError
            raise NCBIQueryError(0, "No orthologs available for NCBI hgid {0}".format(hgid))
        else:
            # If no exception raised, some orthologs exist, raise different NCBIQueryError
            error_msg = "No AGS ortholog present for provided Human Gene ID: {0}".format(hgid)
            raise NCBIQueryError(1, error_msg)
    else:
        spec_href = spec_elem.get_attribute('href')
        spec_gid = re.search("/gene/(\d*)", spec_href).groups()[0]
        id_df.loc[idx,out_col] = spec_gid
    if id_out_fpath:
        id_df.to_csv(id_out_fpath,sep='\t')

def map_AGS_geneIDs(id_df, results_outpath, errors_fpath, tax_dict,
                    overwrite_gid=[]):
    """Reads a csv file at specified path, generates a new DataFrame containing AGS Gene IDs when available.

    :param id_df: DataFrame with gene symbol and NCBI Gene ID information.
     REQUIRED column human_gene_id, empty values okay (but will not fetch ortholog gene IDs for those entries)
    For each row in the tsv file, attempts to fetch NCBI Gene ID for AGS using ortholog information. Will
    fail if 1) Human Gene ID is missing or 2) No orthologs could be fetched for AGS and will write a brief
    error message to errors_tsv.
    :param results_outpath: file path to DataFrame where mapped ortholog Gene IDs will be stored.
    Returns a DataFrame object corresponding to the original csv file with all possible entries of FIELD_NAME
    populated with the ortholog Gene ID numbers.
    :param errors_fpath:
    :param tax_dict:
    :param overwrite_gid:
    :return:
    """
    driver = headless_driver()
    gene_field_name = tax_dict['gid_column']
    #Name of column to store to in DataFrame. Doesn't need to be present in infile.
    if os.path.exists(results_outpath):
        #If results df exists, will default to reading fields/ determining missing entries from there
        out_id_df = pd.read_csv(results_outpath, dtype=str,index_col="overall_index",sep='\t')
    else:
        out_id_df = id_df.copy()
        out_id_df.insert(loc=len(out_id_df.columns),column=gene_field_name)
    check_error_file, NCBI_errors_df = load_errors(errors_fpath, error_type="NCBIQueryError")
    #Determine rows missing gene_field_name and rows missing human_gene_id
    # if gene_field_name not in out_id_df.columns:
    #     missing_spec_gid = out_id_df
    # else:
    missing_spec_gid = out_id_df.loc[out_id_df[gene_field_name].isnull(), :]
    missing_hgid = missing_spec_gid.loc[missing_spec_gid["human_gene_id"].isnull(), :]
    #Iterate over id DataFrame, query NCBI Gene Orthologs if overwrite_gid or if missing
    for idx, row in out_id_df.iterrows():
        symbol = row["gene_symbol"]
        hgid = row["human_gene_id"]
        #Query if force overwrite of data using overwrite_gid or if missing
        if symbol in overwrite_gid or idx in missing_spec_gid.index:
            if check_error_file and symbol in NCBI_errors_df["gene_symbol"].unique():
                print_errors(NCBI_errors_df,symbol)
            elif idx in missing_hgid.index:
                rd_error = RecordDataError(1,"No Human GeneID present in data")
                write_errors(errors_fpath,symbol,rd_error)
                continue
            else:
                try:
                    single_NCBI_gid_query(driver,out_id_df,idx,hgid,tax_dict,id_out_fpath=results_outpath)
                except NCBIQueryError as ncbi_error:
                    write_errors(errors_fpath,symbol,ncbi_error)
    return out_id_df


def download_NCBI_records(id_df, NCBI_records_dirpath,tax_dict, NCBI_API_key,
                          overwrite_fasta=[],pid_outpath=''):
    """Downloads NCBI protein records for each NCBI Gene ID listed in ags_mapped_id_df.

    :param: id_df: DataFrame object with required columns 'Gene Symbol' and 'AGS Gene ID.' Gene symbol
    entries are used to name fasta files downloaded; AGS Gene IDs are queried using Entrez elink
    to 1) match Gene ID to all corresponding Protein IDs and 2) download those Protein IDs into one fasta
    file per gene symbol, saved into NCBI_records_dirpath
    :param NCBI_records_dirpath: directory path where NCBI input data will be downloaded
    :param gene_field_name: column name in ags_mapped_id_df where ortholog IDs are stored
    :param protein_field_name: column name in ags_mapped_id_df where linked Protein ID(s) will be stored as comma separated list
    :param NCBI_API_key: Must be valid API key or empty string

    :return: modified ags_mapped_id_df containing a populated protein_field_name column with linked Protein IDs. Will download
    corresponding Protein Sequences to files named by gene symbol into directory specified by NCBI_records_dirpath
    """
    import subprocess
    import requests
    import xml.etree.ElementTree as ET

    ENTREZ_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    api_key_ext = ""
    if NCBI_API_key:
        #Default value will be empty string if NCBIAPIKey field is empty or absent in config file, defaults
        #to making API requests without API key
        api_key_ext = "&api_key={0}".format(NCBI_API_key)
    gene_field_name, protein_field_name = tax_dict['gid_column'],tax_dict['pid_column']
    ags_mapped_df = id_df.loc[~id_df[gene_field_name].isnull(), :]
    # Convert Gene ID to list of Protein IDs corresponding to transcript variant sequences
    for idx, row in ags_mapped_df.iterrows():
        symbol = row["gene_symbol"]
        fasta_fpath = "{0}/{1}.fasta".format(NCBI_records_dirpath, symbol)
        if not os.path.exists(fasta_fpath) or symbol in overwrite_fasta:
            AGS_gid = row[gene_field_name]
            # print("==={0}===".format(symbol))
            # print("AGS Gene ID: {0}".format(AGS_gid))
            elink_req = "elink.fcgi?dbfrom=gene&db=protein&id={0}{1}".format(AGS_gid, api_key_ext)
            gp_elink_url = ENTREZ_BASE_URL + elink_req

            elink_response = requests.get(gp_elink_url)
            xml_data = elink_response.content
            root = ET.fromstring(xml_data)
            # Check XML formatting of elink pages - update xpath accordingly if functionality breaks
            # Pulls Record IDs for Protein specifically; use gene_protein_refseq for Protein RefSeqs
            protein_IDs = [link.text for link in root.findall(".//LinkSetDb[LinkName='gene_protein']/Link/Id")]
            id_str = ','.join(protein_IDs)
            id_df.loc[idx, protein_field_name] = id_str

            efetch_req = "efetch.fcgi?db=protein&id={0}&rettype=fasta&retmode=text{1}".format(id_str,api_key_ext)
            efetch_url = ENTREZ_BASE_URL + efetch_req
            subprocess.run(args=['wget', efetch_url, '-O', fasta_fpath])
            if pid_outpath:
                id_df.to_csv(pid_outpath, sep='\t')
    return id_df


