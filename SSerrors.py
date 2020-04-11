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