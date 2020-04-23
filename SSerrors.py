# Error classes for OrthoDB query failures, missing sequence data, or GeneCards queries
import os
import pandas as pd
class Error(Exception):
    #Base class
    pass
class SequenceDataError(Error):
    """Error class to be raised if missing GS sequences from OrthoDB/ NCBI data"""
    error_type = "SequenceDataError"
    def __init__(self, code, message):
        self.code = code
        self.message = message


class OrthoDBQueryError(Error):
    """Error class to be raised if OrthoDB query failed to generate input files"""
    error_type = "OrthoDBQueryError"

    def __init__(self, code, message):
        self.code = code
        self.message = message

class RecordDataError(Error):
    """Error class to be raised if NCBI query failed to generate input files
    Codes:
    0 - Missing File
    1 -  Missing human Gene ID from config["IDFilePath"]
    """
    error_type = "RecordDataError"

    def __init__(self, code, message):
        self.code = code
        self.message = message

class NCBIQueryError(Error):
    """Error class to be raised if NCBI query failed to generate input files"""
    error_type = "NCBIQueryError"
    def __init__(self, code, message):
        self.code = code
        self.message = message

class GeneCardsError(Error):
    """Error class if aliases for a gene symbol could not be automatically fetched"""
    error_type = "GeneCardsError"

    def __init__(self, code, message):
        self.code = code
        self.message = message

class SequenceAnalysisError(Error):
    """Error class if JSD/ BLOSUM metrics analysis cannot be completed for a gene"""
    error_type = "SequenceAnalysisError"

    def __init__(self, code, message):
        self.code = code
        self.message = message


def write_errors(errors_fpath, gene_name, error):
    """Maintains a tsv file of gene symbols and errors generated during the run.
    """
    error_type = error.error_type
    error_code = error.code
    error_msg = error.message
    if not os.path.exists(errors_fpath):
        errors_f = open(errors_fpath, 'wt')
        errors_f.write("gene\terror_type\terror_code\terror_str\n")
    else:
        errors_df = pd.read_csv(errors_fpath, delimiter='\t')
        if gene_name in errors_df["gene"].unique():
            gene_error_df = errors_df.loc[errors_df["gene"] == gene_name, :]
            if gene_error_df["error_str"].str.contains(error_msg).any():
                #Checks if error has been logged already, prints to output instead of writing to file
                error_row = gene_error_df.loc[gene_error_df["error_str"] == error_msg, :]
                gene_name, error_type, error_code, error_msg = error_row.values[0]
                print("{0}\t{1}\t{2}\t{3}".format(gene_name, error_type, error_code, error_msg))
                return
    errors_f = open(errors_fpath, 'at')
    fline = "{0}\t{1}\t{2}\t{3}\n".format(gene_name, error_type, error_code, error_msg)
    errors_f.write(fline)
    print(fline)
    errors_f.close()

def print_errors(errors_df,gene_symbol):
    error_row = errors_df.loc[errors_df["gene"] == gene_symbol, :]
    genename, error_type, error_code, error_msg = error_row.values[0]
    print("{0}\t{1}\t{2}\t{3}".format(genename, error_type, error_code, error_msg))

def load_errors(errors_fpath,error_type=""):
    """Loads errors_df from specified file path. Used to prevent repeat operations that generate errors (ie failed
    OrthoDB or NCBI queries, lack of available GeneCards alias data

    :param errors_fpath: file path to errors.tsv which logs errors raised in acquisition/ analysis process. If path
    doesn't exist, returns check_error_file as False and an empty DataFrame (which won't be accessed).
    :param error_type: if provided, will only check against logged errors of given type
    :return: check_error_file: boolean on whether errors_df exists
    :return: errors_df: if check_error, returns file loaded from
    """
    if os.path.exists(errors_fpath):
        errors_df = pd.read_csv(errors_fpath, delimiter='\t')
        if error_type:
            errors_df = errors_df.loc[errors_df["error_type"] == error_type, :]
        check_error_file = True
        return check_error_file, errors_df
    else:
        check_error_file = False
        return check_error_file, pd.DataFrame()

