import unittest
import pandas as pd
import os
import warnings
os.chdir("..")

from SSutility import SSdirectory, SSconfig, SSfasta, SSerrors
from SSacquisition import ODBquery,NCBIquery

test_tmp_dir = "test/tmp"


class testODBQuery(unittest.TestCase):

    def setUp(self):
        SSdirectory.create_directory("{0}/input/ODB".format(test_tmp_dir))

    def test_ODB_Query(self):
        test_spec_str = "species=9606,43179,9601,10090,10161"
        test_level_str = "level=40674"
        #Invalid ODB query
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query(test_tmp_dir,"jsakdh",test_level_str,test_spec_str)
        #ODB Query valid gene symbol no results
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query(test_tmp_dir,"CBR3-AS1",test_level_str,test_spec_str)
        #ODB query valid gene symbol too many clusters
        with self.assertRaises(SSerrors.OrthoDBQueryError):
            ODBquery.ODB_query(test_tmp_dir,"CDC42",test_level_str,test_spec_str)
        #Valid ODB query
        ODBquery.ODB_query(test_tmp_dir, "ATP5MC1", test_level_str, test_spec_str)

        cbr_tsv_path = "{0}/input/ODB/CBR3-AS1.tsv".format(test_tmp_dir)
        cdc42_tsv_path = "{0}/input/ODB/CDC42.tsv".format(test_tmp_dir)
        atp5mc1_tsv_path = "{0}/input/ODB/ATP5MC1.tsv".format(test_tmp_dir)
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
        from IPython.display import display
        # config, spec_list, tax_subset, gene_id_df, tax_table = SSconfig.config_initialization()
        gene_id_df = SSconfig.read_geneID_file("test/test_data/cDNAscreen_geneIDs_clean.csv")
        with pd.option_context('display.max_columns', None, 'display.max_rows', 5):
            pass
        row_generator = gene_id_df.iterrows()
        idx, row = row_generator.__next__()
        symbol = row["gene_symbol"]
        hgid = row["human_gene_id"]
        driver = NCBIquery.headless_driver()
        try:
            test_outfpath = "{0}/mapped_gids.tsv".format(test_tmp_dir)
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

if __name__ == '__main__':

    from SSutility.SSdirectory import create_directory,empty_directory
    create_directory(test_tmp_dir)
    empty_directory(test_tmp_dir)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=[ImportWarning,DeprecationWarning])
        # unittest.main(buffer=False)