import os
import pandas as pd
from SSerrors import NCBIQueryError, write_errors
import re

def download_AGS_data(input_csv_fpath, config, errors_tsv_fpath):
    """Uses NCBI Rest API to acquire sequence data for NCBI gene IDs listed at input_csv_path

    :param input_csv_fpath: file path to table with required column human_gene_id
    :param config: configparser object created from config/config.txt
    :param errors_tsv_fpath: filepath to errors table (see SSerrors.write_errors for details)
    :return: ags_geneID_df: DataFrame containing all columns originally in input_csv_fpath as well as columns
    specified in the config file by NCBIGeneIDField and NCBIProteinIDField with values populated by map_AGS_geneIDs
    and download_NCBI_records respectively.
    """
    from SSdirectory import create_directory
    run_name = config["RunName"]
    NCBI_taxid = config["NCBITaxID"]
    gene_field_name = config["NCBIGeneIDField"]
    protein_field_name = config["NCBIProteinIDField"]
    NCBI_API_key = config.get("NCBIAPIKey","")
    NCBI_input_dir = "{0}/input/NCBI/{1}".format(run_name,NCBI_taxid)
    create_directory(NCBI_input_dir)

    filled_outpath = "{0}/summary/cDNAscreen_geneIDs_complete.csv".format(run_name)
    NCBI_errors_fpath = errors_tsv_fpath

    geneID_df = map_AGS_geneIDs(input_csv_fpath, filled_outpath, NCBI_errors_fpath,gene_field_name,NCBI_taxid)
    ags_geneID_df = geneID_df.loc[~geneID_df[gene_field_name].isnull(), :]
    ags_geneID_df = download_NCBI_records(ags_geneID_df, NCBI_input_dir,gene_field_name,protein_field_name,NCBI_API_key)
    return ags_geneID_df

def NCBI_ID_query(driver,idx,id_df,hgid, query_tax_id="9999",field_name=""):
    """For a human NCBI gene id, searches for an ortholog gene id corresponding to the query_tax_id

    :param driver: Selenium webdriver object, initialized pre-iteration with headless option
    :param idx: for row in id_df in which to store the ortholog gene ID
    :param id_df: DataFrame object where ortholog gene IDs will be stored
    :param query_tax_id: NCBI Taxonomy ID corresponding to species for which orthologs will be fetched.
    :param field_name: Column in id_df in which orthologs will be stored. Will be created if not present. Defaults to
            [query_tax_id]_gene_ID if no value is provided.
    """
    from selenium.common.exceptions import NoSuchElementException
    if not field_name:
        field_name = "{0}_gene_id".format(query_tax_id)
    try:
        req_url = "https://www.ncbi.nlm.nih.gov/gene/{0}/ortholog/?scope={1}".format(hgid,query_tax_id)
        driver.get(req_url)
        result_url = driver.current_url
        if re.search("scope={0}".format(query_tax_id), result_url):
            entry_xpath = "//tbody/tr/td[@class=' fld-gene']/a"
            entry = driver.find_element_by_xpath(entry_xpath)
            entry_href = entry.get_attribute('href')
            entry_gid = re.search("/gene/(\d*)", entry_href).groups()[0]
            id_df.loc[idx, field_name] = entry_gid
        else:
            error_msg = "No AGS ortholog present for provided Human Gene ID: {0}".format(hgid)
            raise NCBIQueryError(0,error_msg)
    except NoSuchElementException as e:
        error_msg = "No Orthologs present for GeneID {0}".format(hgid)
        raise NCBIQueryError(1,error_msg)


def map_AGS_geneIDs(csv_inpath, results_outpath, errors_tsv_path, gene_field_name, query_tax_id):
    """Reads a csv file at specified path, generates a new DataFrame containing AGS Gene IDs when available.

    csv_inpath: csv file containing list of genes and other identifying information. REQUIRED column human_gene_id,
    empty values okay (but will not fetch ortholog gene IDs for those entries)
    For each row in the csv file, attempts to fetch NCBI Gene ID for AGS using ortholog information. Will
    fail if 1) Human Gene ID is missing or 2) No orthologs could be fetched for AGS and will write a brief
    error message to errors_tsv.
    Returns a DataFrame object corresponding to the original csv file with all possible entries of FIELD_NAME
    populated with the ortholog Gene ID numbers.

    """

    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    WINDOW_SIZE = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % WINDOW_SIZE)
    driver = webdriver.Chrome(chrome_options=chrome_options)

    #Name of column to store to in DataFrame. Doesn't need to be present in infile.
    field_conv_dict = {"human_gene_id": str, gene_field_name: str}
    if os.path.exists(results_outpath):
        #If results df exists, will default to reading fields/ determining missing entries from there
        id_df = pd.read_csv(results_outpath, dtype=field_conv_dict)
    else:
        id_df = pd.read_csv(csv_inpath, dtype=field_conv_dict)
    #Determine rows missing gene_field_name and rows missing human_gene_id
    missing_AGS_gid = id_df.loc[id_df[gene_field_name].isnull(), :]
    missing_hgid = missing_AGS_gid.loc[missing_AGS_gid["human_gene_id"].isnull(), :]
    #For any row missing gene_field_name, will
    for idx, row in missing_AGS_gid.iterrows():
        symbol = row["gene_symbol"]
        hgid = row["human_gene_id"]
        if idx in missing_hgid.index:
            write_errors(symbol, "No Human GeneID present in data", errors_tsv_path)
            continue
        else:
            try:
                NCBI_ID_query(driver,idx,hgid,query_tax_id,field_name=gene_field_name)
            except NCBIQueryError as ncbi_error:
                write_errors(errors_tsv_path,symbol,ncbi_error)

    id_df.to_csv(results_outpath)
    return id_df


def download_NCBI_records(ags_geneID_df, NCBI_records_dirpath,gene_field_name, protein_field_name, NCBI_API_key):
    """Downloads NCBI protein records for each NCBI Gene ID listed in ags_geneID_df.

    :param: ags_geneID_df: DataFrame object with required columns 'Gene Symbol' and 'AGS Gene ID.' Gene symbol
    entries are used to name fasta files downloaded; AGS Gene IDs are queried using Entrez elink
    to 1) match Gene ID to all corresponding Protein IDs and 2) download those Protein IDs into one fasta
    file per gene symbol, saved into NCBI_records_dirpath
    :param NCBI_records_dirpath: directory path where NCBI input data will be downloaded
    :param gene_field_name: column name in ags_geneID_df where ortholog IDs are stored
    :param protein_field_name: column name in ags_geneID_df where linked Protein ID(s) will be stored as comma separated list
    :param NCBI_API_key: Must be valid API key or empty string

    :return: modified ags_geneID_df containing a populated protein_field_name column with linked Protein IDs. Will download
    corresponding Protein Sequences to files named by gene symbol into directory specified by NCBI_records_dirpath
    """
    import subprocess
    import urllib.request
    import xml.etree.ElementTree as ET

    ENTREZ_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    api_key_url_ext = ""
    #You can add gene symbols from your gene list to OVERWERITE_FASTAS if you want the data to be redownloaded/ updated
    #every time the script is run
    OVERWRITE_FASTAS = []
    if NCBI_API_key:
        #Default value will be empty string if NCBIAPIKey field is empty or absent in config file, defaults
        #to making API requests without API key
        api_key_url_ext = "&api_key={0}".format(NCBI_API_key)
    # Convert Gene ID to list of Protein IDs corresponding to transcript variant sequences
    for idx, row in ags_geneID_df.iterrows():
        symbol = row["Gene Symbol"]
        print("==={0}===".format(symbol))
        fasta_fpath = "{0}/{1}_AGS.fasta".format(NCBI_records_dirpath, symbol)
        if not os.path.exists(fasta_fpath) or symbol in OVERWRITE_FASTAS:
            AGS_gid = row[gene_field_name]
            print("AGS Gene ID: {0}".format(AGS_gid))
            elink_req = "elink.fcgi?dbfrom=gene&db=protein&id={0}{1}".format(AGS_gid, api_key_url_ext)
            gp_elink_url = ENTREZ_BASE_URL + elink_req

            file = urllib.request.urlopen(gp_elink_url)
            xml_data = file.read()
            file.close()

            root = ET.fromstring(xml_data)
            # Check XML formatting of elink pages - update xpath accordingly if functionality breaks
            # Pulls Record IDs for Protein specifically; use gene_protein_refseq for Protein RefSeqs
            protein_IDs = [link.text for link in root.findall(".//LinkSetDb[LinkName='gene_protein']/Link/Id")]
            id_str = ','.join(protein_IDs)
            ags_geneID_df.loc[idx, protein_field_name] = id_str

            efetch_req = "efetch.fcgi?db=protein&id={0}&rettype=fasta&retmode=text{1}".format(id_str,api_key_url_ext)
            efetch_url = ENTREZ_BASE_URL + efetch_req
            subprocess.run(args=['wget', efetch_url, '-O', fasta_fpath])
    return ags_geneID_df