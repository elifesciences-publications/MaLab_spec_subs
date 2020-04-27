import SSconfig
import SSdirectory
import SSerrors


def ss_acquisition(config, spec_list, tax_subet, gene_id_df, tax_table):
    import ODBquery
    import NCBIquery
    import aliasQuery

    gene_symbols = gene_id_df['gene_symbol']
    valid_queries, failed_queries = ODBquery.download_ODB_input(gene_symbols, tax_table, config)
    ags_mapped_df = NCBIquery.download_AGS_data(gene_id_df, config)
    aliasQuery.download_alias_data(gene_symbols, config)

def ss_filter(config, spec_list, tax_subet, gene_id_df, tax_table):



def main():
    config, spec_list, tax_subet, gene_id_df, tax_table = SSconfig.config_initialization()
    ss_acquisition(config, spec_list, tax_subet, gene_id_df, tax_table)
    ss_filter(config, spec_list, tax_subet, gene_id_df, tax_table)
