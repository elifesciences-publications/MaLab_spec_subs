import os
import pandas as pd
from SSerrors import write_errors, print_errors, load_errors, OrthoDBQueryError
# Acquire input data via OrthoDB API
def ODB_query(run_name, gene_name, level_str, spec_str):
    """Queries OrthoDB via the fasta and tab API for gene_name.
    More info: https://www.orthodb.org/orthodb_userguide.html#api
    level_str corresponds to the API variable for phylogenetic clade
    spec_str corresponds to the taxonomy ids for the list of species from the config folder
    """
    import time, json
    from json import JSONDecodeError
    from bs4 import BeautifulSoup
    import subprocess
    # File paths and OrthoDB urls for downloads. NOTE BASE_URL might need updating depending on ODB conventions
    BASE_URL = "https://v101.orthodb.org"
    query_str = "query={0}".format(gene_name)
    fasta_url = "{0}/fasta?{1}&{2}&{3}".format(BASE_URL, query_str, level_str, spec_str)
    fasta_path = "{0}/input/ODB/{1}.fasta".format(run_name, gene_name)
    tsv_url = "{0}/tab?{1}&{2}&{3}".format(BASE_URL, query_str, level_str, spec_str)
    tsv_path = "{0}/input/ODB/{1}.tsv".format(run_name, gene_name)
    # Obey OrthoDB download restrictions (one request per second) bc you're a good noodle
    t1 = time.process_time()
    fasta_proc = subprocess.run(args=['wget', fasta_url, '-O', fasta_path])
    if (time.process_time() - t1) < 1:
        time.sleep(0.5)
    t1 = time.process_time()
    tsv_proc = subprocess.run(args=['wget', tsv_url, '-O', tsv_path])
    if (time.process_time() - t1) < 1:
        time.sleep(0.5)
    try:
        # JSON format returned if no results for query string - try opening downloaded data as JSON, if
        # successful, raise an OrthoDBQueryError
        tsv_json = json.load(open(tsv_path))
        os.remove(fasta_path)
        os.remove(tsv_path)
        raise OrthoDBQueryError(0, "No OrthoDB results for query")
    except JSONDecodeError:
        # Check if html syntax present in file (result of too many clusters returned to be downloaded);
        # if not, query was successful and run_name/input should now have ODB formatted .fasta and .tsv files
        file_txt = ""
        with open(fasta_path, "rt") as fasta_f:
            for i in range(10):
                file_txt = file_txt + fasta_f.readline()
            if bool(BeautifulSoup(file_txt, "html.parser").find()):
                os.remove(fasta_path)
                os.remove(tsv_path)
                raise OrthoDBQueryError(1, "OrthoDB search yielded too many clusters")
                # If no OrthoDBQueryError is raised, download was successful (no further action needed)


def download_ODB_input(gene_list, tax_table, config):
    """Queries OrthoDB for all entries in gene list (logs failed searches into errors_fpath), using species
    list from tax_table and taxonomy level provided in config directory. This function will attempt to
    query OrthoDB for each gene symbol in gene_list according to the species list in the config directory.
    Note that for configuring the species list, one OrthoDB input file can be generated with a large
    number of species whose sequences can later be removed from the sequence set used for alignment/ analysis,
    allowing the same input file to serve different downstream analyses with smaller species sets.

    :param gene_list: iterable (ie list or Pandas Series) containing gene symbols for which OrthoDB data
    will be downloaded
    :param tax_table: Table constructed from OrthoDB raw species table, contains species name, OrthoDB ID, and
    assembly information.
    :param config: configparser object constructed from config/config.txt, contains run parameters

    Returns the list of gene symbols from gene_list for which OrthoDB data was successfully downloaded
    and the list of gene symbols for which the OrthoDB queries failed"""
    tax_ids = tax_table["tax_id"].values.astype(str)
    spec_str = "species=" + ",".join(tax_ids)
    level_str = "level=" + str(config["ODBLevel"])
    failed_queries = []
    run_name = config["RunName"]
    errors_fpath = config["ErrorsFilePath"]
    check_error_file, ODB_errors_df = load_errors(errors_fpath,error_type="OrthoDBQueryError")
    # if os.path.exists(errors_fpath):
    #     errors_df = pd.read_csv(errors_fpath, delimiter='\t')
    #     ODB_errors_df = errors_df.loc[errors_df["error_type"] == "OrthoDBQueryError", :]
    #     check_error_file = True
    # else:
    #     check_error_file = False
    for gene_name in gene_list:
        fasta_path = "{0}/input/ODB/{1}.fasta".format(run_name, gene_name)
        if config.getboolean("OverwriteInput") or not os.path.exists(fasta_path):
            if check_error_file and gene_name in ODB_errors_df["gene"].unique():
                print_errors(ODB_errors_df,gene_name)
                failed_queries.append(gene_name)
            else:
                try:
                    ODB_query(run_name, gene_name, level_str, spec_str)
                except OrthoDBQueryError as odb_error:
                    failed_queries.append(gene_name)
                    write_errors(errors_fpath, gene_name, odb_error)

    print("Input queries downloaded.")
    valid_queries = [gene for gene in gene_list if gene not in failed_queries]
    return valid_queries, failed_queries
