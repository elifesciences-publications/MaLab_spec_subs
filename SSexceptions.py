# Exception classes for OrthoDB query failures, missing sequence data, or GeneCards queries
class Exception(Exception):
    """Base class for exceptions in this module."""
    pass


class SequenceDataException(Exception):
    """Exception class to be raised if missing GS sequences from OrthoDB/ NCBI data"""
    exception_type = "SequenceDataException"

    def __init__(self, code, message):
        self.code = code
        self.message = message


class OrthoDBQueryException(Exception):
    """Exception class to be raised if OrthoDB query failed to generate input files"""
    exception_type = "OrthoDBQueryException"

    def __init__(self, code, message):
        self.code = code
        self.message = message


class GeneCardsException(Exception):
    """Exception class if aliases for a gene symbol could not be automatically fetched"""
    exception_type = "GeneCardsException"

    def __init__(self, code, message):
        self.code = code
        self.message = message


class SequenceAnalysisException(Exception):
    """Exception class if JSD/ BLOSUM metrics analysis cannot be completed for a gene"""
    exception_type = "SequenceAnalysisException"

    def __init__(self, code, message):
        self.code = code
        self.message = message


def write_exceptions(exceptions_fpath, gene_name, exception):
    """Maintains a tsv file of gene symbols and exceptions generated during the run.
    """
    exception_type = exception.exception_type
    exception_code = exception.code
    exception_msg = exception.message
    if not os.path.exists(exceptions_fpath):
        exceptions_f = open(exceptions_fpath, 'wt')
        exceptions_f.write("gene\texception_type\texception_code\texception_str\n")
    else:
        exceptions_df = pd.read_csv(exceptions_fpath, delimiter='\t')
        if gene_name in exceptions_df["gene"].unique():
            gene_exception_df = exceptions_df.loc[exceptions_df["gene"] == gene_name, :]
            if gene_exception_df["exception_str"].str.contains(exception_msg).any():
                #                 print("Previously stored exception:")
                exception_row = gene_exception_df.loc[gene_exception_df["exception_str"] == exception_msg, :]
                gene_name, exception_type, exception_code, exception_msg = exception_row.values[0]
                print("{0}\t{1}\t{2}\t{3}".format(gene_name, exception_type, exception_code, exception_msg))
                return
    exceptions_f = open(exceptions_fpath, 'at')
    fline = "{0}\t{1}\t{2}\t{3}\n".format(gene_name, exception_type, exception_code, exception_msg)
    exceptions_f.write(fline)
    print(fline)
    exceptions_f.close()