import unittest
import os
os.chdir("..")
from SSutility import SSdirectory, SSerrors

import contextlib
import io
import pandas as pd
from SSfilter import ODBfilter

test_temp_dir = "test/tmp"

class SSErrorsTest(unittest.TestCase):

    def setUp(self):
        SSdirectory.empty_directory(test_temp_dir)

    def test_write_errors(self):
        test_efpath = "{0}/test_errors_fpath.tsv".format(test_temp_dir)
        self.assertFalse(os.path.exists(test_efpath))
        print("Inserting test error")
        test_error = SSerrors.NCBIQueryError(0,"test_message")
        SSerrors.write_errors(test_efpath,"CD151",test_error)
        self.assertTrue(os.path.exists(test_efpath))
        check_ef, errors_df = SSerrors.load_errors(test_efpath)
        self.assertTrue(len(errors_df) == 1)
        print("reinserting same error")
        SSerrors.write_errors(test_efpath, "CD151", test_error)
        check_ef, errors_df = SSerrors.load_errors(test_efpath)
        self.assertTrue(len(errors_df) == 1)
        print("Adding second error")
        test_error2 = SSerrors.NCBIQueryError(0,"different test message")
        SSerrors.write_errors(test_efpath, "CD151", test_error2)
        check_ef, errors_df = SSerrors.load_errors(test_efpath)
        self.assertTrue(len(errors_df)==2)


    def test_print_errors(self):
        test_efpath = "{0}/test_errors_fpath.tsv".format(test_temp_dir)

        test_error = SSerrors.NCBIQueryError(0, "test_message")
        SSerrors.write_errors(test_efpath, "CD151", test_error)
        test_error2 = SSerrors.NCBIQueryError(0, "different test message")
        SSerrors.write_errors(test_efpath, "CD151", test_error2)

        check,errors_df = SSerrors.load_errors(test_efpath)
        out_buf = io.StringIO()
        with contextlib.redirect_stdout(out_buf):
            SSerrors.print_errors(errors_df,"CD151")
            output = out_buf.getvalue()
            for i,row in errors_df.iterrows():
                genename, error_type, error_code, error_msg = row.values
                out_line = ("{0}\t{1}\t{2}\t{3}".format(genename, error_type, error_code, error_msg))
                self.assertIn(out_line,output)
        #Test message specific printing
        out_buf = io.StringIO()
        with contextlib.redirect_stdout(out_buf):
            SSerrors.print_errors(errors_df, "CD151","test_message")
            output = out_buf.getvalue()
            for i, row in errors_df.iterrows():
                genename, error_type, error_code, error_msg = row.values
                out_line = ("{0}\t{1}\t{2}\t{3}".format(genename, error_type, error_code, error_msg))
                if i == 0:
                    self.assertTrue(out_line in output)
                else:
                    self.assertFalse(out_line in output)

    def test_qc_log(self):
        test_fpath = "{0}/qc_test.tsv".format(test_temp_dir)
        #Write one entry
        ODBfilter.write_ref_seq_QC(test_fpath,"CD151","a test message")
        test_qc = pd.read_csv(test_fpath,sep="\t",index_col=0)
        self.assertTrue(len(test_qc) == 1)
        #Attempt rewriting same entry
        ODBfilter.write_ref_seq_QC(test_fpath, "CD151", "a test message")
        test_qc = pd.read_csv(test_fpath, sep="\t", index_col=0)
        self.assertTrue(len(test_qc) == 1)
        # Test new entry writing with same symbol identifier
        ODBfilter.write_ref_seq_QC(test_fpath, "CD151", "a new test message!")
        test_qc = pd.read_csv(test_fpath, sep="\t", index_col=0)
        self.assertTrue(len(test_qc) == 2)


    def test_load_errors(self):
        #Checks file path for
        from SSutility import config
        from IPython.display import display
        errors_fpath = config['RUN']['ErrorsFilePath']
        print(errors_fpath)
        check, errors_df = SSerrors.load_errors(errors_fpath)
        display(errors_df)
        self.assertTrue("gene_symbol" in errors_df.columns)

if __name__ == "__main__":
    unittest.main()