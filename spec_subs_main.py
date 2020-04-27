from SSutility import SSconfig, SSdirectory, SSerrors
import sys, os

def ss_acquisition(config, spec_list, tax_subet, gene_id_df, tax_table):
    # import ODBquery
    # import NCBIquery
    # import aliasQuery
    from SSacquisition import ODBquery,NCBIquery,aliasQuery
    gene_symbols = gene_id_df['gene_symbol']
    valid_queries, failed_queries = ODBquery.download_ODB_input(gene_symbols, tax_table, config)
    ags_mapped_df = NCBIquery.download_AGS_data(gene_id_df, config)
    aliasQuery.download_alias_data(gene_symbols, config)

def ss_filter(config, spec_list, tax_subet, gene_id_df, tax_table):

    from SSfilter import ODBfilter,NCBIfilter
    run_name,errors_fpath= config['RUN']['RunName'],config['RUN']['ErrorsFilePath']
    ncbi_taxid,ncbi_specname = config['NCBI']['NCBITaxID'],config['NCBI']['NCBITaxName']
    odb_test_taxid = config['ODB']['ODBTestTaxID']
    gene_symbols = gene_id_df['gene_symbol']
    tax_subset = list(config['AnalysisODBTaxSubset'].keys())
    taxid_dict = {ncbi_specname: ncbi_taxid}

    check_errors,errors_df = SSerrors.load_errors(errors_fpath)
    for symbol in gene_symbols:
        out_records_fpath = "{0}/output/{1}/{1}_records.tsv".format(run_name, symbol)
        if not os.path.exists(out_records_fpath) or config['RUN'].getboolean('OverwriteFilter'):
            if check_errors and symbol in errors_df['gene_symbol'].unique():
                odb_error = (symbol in errors_df.loc[errors_df['error_type']=="OrthoDBQueryError",'gene_symbol'].unique())
                ncbi_error = (symbol in errors_df.loc[errors_df['error_type'] == "NCBIQueryError", 'gene_symbol'].unique())
                if odb_error or ncbi_error:
                    continue
            try:
                odb_fpath = "{0}/input/ODB/{1}.fasta".format(run_name,symbol)
                ncbi_fpath = "{0}/input/NCBI/{1}/{2}.fasta".format(run_name,ncbi_taxid,symbol)
                results = ODBfilter.process_input(symbol, config, tax_subset)
                final_odb, em_df, am_df = results['final_df'], results['em_df'], results['am_df']
                final_combined = NCBIfilter.select_NCBI_record(odb_fpath, ncbi_fpath, taxid_dict,
                                                               final_odb, [odb_test_taxid])
                comb_proc = NCBIfilter.combined_records_processing(config, am_df, em_df, final_combined, symbol)
            except SSerrors.SequenceDataError as e:
                continue

def main():
    import SSutility
    from SSutility import config, spec_list, tax_subet, gene_id_df, tax_table
    # ss_acquisition(config, spec_list, tax_subet, gene_id_df, tax_table)
    ss_filter(config, spec_list, tax_subet, gene_id_df, tax_table)

if __name__ == '__main__':
    main()
