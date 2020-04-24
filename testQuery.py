import unittest
import pandas as pd
import os

import SSdirectory, SSconfig, SSfasta, SSerrors
import NCBIquery

class testNCBIQuery(unittest.TestCase):


    def test_NCBI_query(self):
        # Ignore me
        import SSconfig
        from IPython.display import display
        config, spec_list, gene_id_df, tax_table = SSconfig.config_initialization()
        with pd.option_context('display.max_columns', None, 'display.max_rows', 5):
            pass
        row_generator = gene_id_df.iterrows()
        idx, row = row_generator.__next__()
        symbol = row["gene_symbol"]
        hgid = row["human_gene_id"]
        driver = NCBIquery.headless_driver()
        try:
            test_outfpath = "tmp/mapped_gids.tsv"
            test_ins_col = '9999_gid'
            tax_dict = {'gid_column':test_ins_col,'spec_name':'Urocitellus parryii'}
            NCBIquery.single_NCBI_gid_query(driver, gene_id_df, idx, hgid,tax_dict,
                                                          id_out_fpath=test_outfpath)
            self.assertTrue(test_ins_col in gene_id_df.columns)
        except SSerrors.NCBIQueryError as ncbi_e:
            print(driver.current_url)
            print(ncbi_e)
            driver.quit()

        driver.quit()