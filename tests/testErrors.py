#testErrors.py - Unit tests for error handling
#and logging helper functions for storing errors.
# Copyright (C) 2020  Evan Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import os
# os.chdir("..")
import sys
sys.path.append(os.getcwd())
from SSutility import SSdirectory, SSerrors

import contextlib
import io
import pandas as pd
from SSfilter import ODBfilter

test_temp_dir = "tests/tmp"

class SSErrorsTest(unittest.TestCase):

    def setUp(self):
        SSdirectory.empty_directory(test_temp_dir)

    def test_write_errors(self):
        print("Write Errors Test")
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
        print("Print Errors Test - adding test errors")
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
        print("Quality check log test")
        test_fpath = "{0}/qc_test.tsv".format(test_temp_dir)
        #Write one entry
        print("Writing test entry, should only output once")
        ODBfilter.write_ref_seq_QC(test_fpath,"CD151","a test message")
        test_qc = pd.read_csv(test_fpath,sep="\t",index_col=0)
        self.assertTrue(len(test_qc) == 1)
        #Attempt rewriting same entry
        print("rewriting same entry - should output only once")
        ODBfilter.write_ref_seq_QC(test_fpath, "CD151", "a test message")
        test_qc = pd.read_csv(test_fpath, sep="\t", index_col=0)
        self.assertTrue(len(test_qc) == 1)
        print("Adding new entry - should only print new entry")
        # Test new entry writing with same symbol identifier
        ODBfilter.write_ref_seq_QC(test_fpath, "CD151", "a new test message.")
        test_qc = pd.read_csv(test_fpath, sep="\t", index_col=0)
        self.assertTrue(len(test_qc) == 2)


    def test_load_errors(self):
        #Checks file path for
        from SSutility import config
        from IPython.display import display
        run_name, errors_fname = config['RUN']['RunName'], config['RUN']['ErrorsFileName']
        errors_fpath = "{0}/{1}".format(run_name, errors_fname)
        print(errors_fpath)
        check, errors_df = SSerrors.load_errors(errors_fpath)
        # display(errors_df)
        self.assertTrue("gene_symbol" in errors_df.columns)

if __name__ == "__main__":
    unittest.main()