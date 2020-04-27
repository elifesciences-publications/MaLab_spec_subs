import unittest
import pandas as pd
import os

import SSdirectory, SSconfig, SSfasta, SSerrors
import ODBquery
import NCBIquery

class testODBQuery(unittest.TestCase):

    def setUp(self):
        SSdirectory.create_directory("tmp/input/ODB")

    def test_ODB_Query(self):
        test_spec_str = "species=9606,43179,9601,10090,10161"
        test_level_str = "level=40674"
        #Invalid ODB query
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query("tmp","jsakdh",test_level_str,test_spec_str)
        #ODB Query valid gene symbol no results
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query("tmp","CBR3-AS1",test_level_str,test_spec_str)
        #ODB query valid gene symbol too many clusters
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query("tmp","CDC42",test_level_str,test_spec_str)
        #Valid ODB query
        ODBquery.ODB_query("tmp", "ATP5MC1", test_level_str, test_spec_str)

        cbr_tsv_path = "tmp/input/ODB/CBR3-AS1.tsv"
        cdc42_tsv_path = "tmp/input/ODB/CDC42.tsv"
        atp5mc1_tsv_path = "tmp/input/ODB/ATP5MC1.tsv"
        #Check that bad ODB queries clear attempted input files
        self.assertFalse(os.path.exists(cbr_tsv_path))
        self.assertFalse(os.path.exists(cdc42_tsv_path))
        self.assertTrue(os.path.exists(atp5mc1_tsv_path))
        #Check downloaded data formatted correctly
        test_tsv = SSfasta.load_tsv_table(atp5mc1_tsv_path)
        self.assertTrue("43179_0" in test_tsv['organism_taxid'].unique())

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