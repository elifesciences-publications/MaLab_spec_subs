# from SSutility import SSconfig, SSdirectory, SSerrors
import sys

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
    pass


def main():
    print(sys.path)
    import SSutility
    from SSutility import config, spec_list, tax_subet, gene_id_df, tax_table
    print(spec_list)
    # SSconfig.config_initialization()
    # ss_acquisition(config, spec_list, tax_subet, gene_id_df, tax_table)
    # ss_filter(config, spec_list, tax_subet, gene_id_df, tax_table)

if __name__ == '__main__':
    main()
