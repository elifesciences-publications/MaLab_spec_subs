#testFilter.py - Unit tests for record selection and data filtering
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

import os
os.chdir("..")

import pandas as pd
import unittest
from SSutility import SSfasta, SSconfig
from IPython.display import display
import warnings
import subprocess
from SSfilter import ODBfilter, NCBIfilter

test_data_dir = "tests/test_data"
test_tmp_dir = "tests/tmp"

 

class SSfastaTest(unittest.TestCase):

    def test_length_load(self):
        test_fpath = "{0}/ODB/{1}.fasta".format(test_data_dir,'ATP5MC1')
        unfiltered_seqs, unfiltered_lens = SSfasta.length_srs(test_fpath)
        self.assertTrue(136 in unfiltered_lens.unique())
        self.assertTrue(138 in unfiltered_lens.unique())
        self.assertTrue('9544_0:0008ab' in unfiltered_lens.index)
        self.assertTrue(len(unfiltered_lens) == 29)
        filtered_ids = ['10090_0:0034c4','43179_0:00103c','9606_0:00415a']
        filtered_seqs, filtered_lens = SSfasta.length_srs(test_fpath,filtered_ids)
        self.assertTrue(136 in filtered_lens.unique())
        self.assertFalse(138 in filtered_lens.unique())
        self.assertFalse('9544_0:0008ab' in filtered_lens.index)
        self.assertTrue(len(filtered_lens) == 3)

    def test_filter_infile(self):
        from Bio import SeqIO
        test_fpath = "{0}/ODB/{1}.fasta".format(test_data_dir,'ATP5MC1')
        ordered_test_ids = ["10090_0:0034c4","43179_0:00103c","9606_0:00415a","10116_0:00386d","42254_0:001ba2",
                            "9986_0:0033f5"]

        unordered_fpath = "{0}/ATP5MC1_unordered.fasta".format(test_tmp_dir)
        ordered_fpath = "{0}/ATP5MC1_ordered.fasta".format(test_tmp_dir)
        SSfasta.filter_fasta_infile(ordered_test_ids,test_fpath,outfile_path=unordered_fpath,ordered=False)
        SSfasta.filter_fasta_infile(ordered_test_ids, test_fpath, outfile_path=ordered_fpath, ordered=True)

        unordered = SeqIO.parse(unordered_fpath, 'fasta')
        ordered = SeqIO.parse(ordered_fpath,'fasta')
        #Check order
        unordered_test_ids = [ordered_test_ids[i] for i in [0, 3, 4, 1, 2, 5]]
        for i,fasta in enumerate(unordered):
            self.assertTrue(fasta.id == unordered_test_ids[i])

        for i,fasta in enumerate(ordered):
            self.assertTrue(fasta.id == ordered_test_ids[i])


class ODBFilterFunctionTest(unittest.TestCase):
    def test_alias_loading(self):
        gene_id_df = SSconfig.read_geneID_file("{0}/cDNAscreen_geneIDs_clean.csv".format(test_data_dir))
        gene_symbols = gene_id_df["gene_symbol"]
        gc_names = []
        for symbol in gene_symbols:
            alias_fpath = "alias_data/{0}_aliases.txt".format(symbol)
            if not os.path.exists(alias_fpath):
                print("Non-existent alias data file path")
                print(alias_fpath)
                continue
            with open(alias_fpath,'r') as alias_f:
                aliases = [alias.strip() for alias in alias_f.readlines()]
                alias_f.close()
            if len(aliases) == 1:
                print("symbol: {0} has only one alias".format(symbol))
            gc_names.append(aliases[0])
        print(len(gc_names))

        self.assertTrue(len(gene_symbols) == 380)
        #Should be missing one entry for CACB1
        self.assertTrue(len(gc_names) == 379)
        self.assertTrue(sum([symbol in gc_names for symbol in gene_symbols]) == 370)
        mismatch_symbols = [symbol for symbol in gene_symbols if symbol not in gc_names]
        self.assertTrue("ATPIF1" in mismatch_symbols and "ISPD" in mismatch_symbols)
    
    
    def test_alias_pattern_match(self):
        from SSfilter.ODBfilter import format_odb_field, odb_field_to_re
        import re
        alias_fpath = "alias_data/ISPD_aliases.txt"
        alias_f = open(alias_fpath,'rt')
        aliases = [alias.strip() for alias in alias_f.readlines()]
        alias_f.close()
        gc_name = aliases[0]
    
        formatted_aliases = [format_odb_field(alias) for alias in aliases]
        alias_res = [odb_field_to_re(alias) for alias in formatted_aliases]
        aliases_pat = "({0})".format("|".join(alias_res))
        test_str_list = ["","CrpPa"," CrPPa   \n","CDP-L-Ribitol Pyrophosphorylase A", "CDP-L-Ribitol Pyrophosphorylase", \
                         "Testicular Tissue Protein Li 97","4-Diphosphocytidyl-2C-Methyl-D-Erythritol Synthase Homolog (Arabidopsis)"]
        test_results = [False,True,True,True,False,True,True]
        for i,test_str in enumerate(test_str_list):
            formatted_test = format_odb_field(test_str)
            result = re.search(aliases_pat, formatted_test) is not None
            try:
                self.assertTrue(result == test_results[i])
            except AssertionError as e:
                print("test str:{0}".format(test_str))
                print("formatted test:{0}".format(formatted_test))
                print("Failed. Result: {0}, Expected result: {1}".format(result, test_results[i]))
                raise e

    def test_alias_df_filter(self):
        from SSfilter.ODBfilter import find_alias_matches
        pre_lengths = [31,29,312,179]
        post_lengths = [31,29,27,125]
        symbol_list = ["ISPD","ATP5MC1","APEX1","CALM1"]
        try:
            for i,symbol in enumerate(symbol_list):
                test_tsv = pd.read_csv("{0}/ODB/{1}.tsv".format(test_data_dir,symbol),sep='\t',index_col='int_prot_id')

                self.assertTrue(len(test_tsv) == pre_lengths[i])
                am_ids, exact_matches = find_alias_matches(symbol,test_tsv,"{0}/errors.tsv".format(test_tmp_dir))

                self.assertTrue(len(am_ids) == post_lengths[i])
                print("Passed {0}/{1}".format(i+1,len(symbol_list)))
        except AssertionError:
            print("Failed for symbol: {0}".format(symbol))
            print("Unfiltered length: {0}, Expected length: {1}".format(len(test_tsv),pre_lengths[i]))
            print("alias match length: {0}, Expected length: {1}".format(len(am_ids),post_lengths[i]))
            print("Alias Matched DF")
            display(test_tsv.loc[am_ids,:])


    def test_kalign(self):
        test_inpath = "{0}/ODB/ATP5MC1.fasta".format(test_data_dir)
        test_outpath = "{0}/test_aln.fasta".format(test_tmp_dir)
        print("Current kalign bin: ")
        subprocess.run(args=["which", "kalign"])
        print("Testing KALIGN output files. ")
        self.assertTrue(os.path.exists(test_inpath))
        from subprocess import Popen, PIPE
        no_args = ["kalign"]
        with open(test_inpath,'r') as test_in, open(test_outpath,'wt',encoding='utf-8') as test_out:
            proc = subprocess.run(args=no_args,stdin=test_in,stdout=test_out,text=True)
        self.assertTrue(os.path.exists(test_outpath))
        with open(test_outpath,'r') as test_out_read:
            self.assertTrue(len(test_out_read.readlines()) > 0)


    def test_exact_match_df(self):
        from SSfilter.ODBfilter import exact_match_df
        test_input_path = "{0}/ODB/ATP5MC1.tsv".format(test_data_dir)
        tsv_df = SSfasta.load_tsv_table(test_input_path)
        unfiltered_uniques = tsv_df['pub_gene_id'].unique()
        self.assertTrue('ATP5G1;ATP5MC1' in unfiltered_uniques)
        self.assertTrue('ATP5G1' in unfiltered_uniques)
        self.assertTrue('ATP5MC1' in unfiltered_uniques)
        self.assertTrue('9823_0:003c30' in tsv_df.index and 'LOC100519871' in unfiltered_uniques)
        test_em = ['ATP5G1','ATP5MC1']
        exact_matches = exact_match_df(tsv_df,test_em)
        pgid_uniques = exact_matches['pub_gene_id'].unique()
        self.assertTrue('ATP5G1;ATP5MC1' in pgid_uniques)
        self.assertTrue('ATP5G1' in pgid_uniques)
        self.assertTrue('ATP5MC1' in pgid_uniques)
        self.assertFalse('LOC100519871' in pgid_uniques)

    # @unittest.skip("broken")
    def test_ksr_record_select(self):
        import SSfilter.ODBfilter
        test_symbol_list = ['ATP5MC1','CALM1','ATPIF1','CD151']
        tax_subset = ['10090_0','43179_0','9606_0','10116_0','42254_0','9601_0']
        errors_fpath =  "cDNAscreen_041020/summary/errors.tsv"
        ks_tids = ['10090_0','43179_0','9606_0']
        # errors_fpath = 'tmp/'
        display_match_data = False
        display_ksr = True

        tmp_manual_selections = "{0}/manual_record_selections.tsv".format(test_tmp_dir)

        with pd.option_context('display.max_columns',None,'display.max_colwidth',500):
            for symbol in test_symbol_list:
                tsv_inpath = "{0}/ODB/{1}.tsv".format(test_data_dir,symbol)
                unfiltered_tsv = SSfasta.load_tsv_table(tsv_inpath,tax_subset=tax_subset)
                unfiltered_fasta = "{0}/ODB/{1}.fasta".format(test_data_dir,symbol)
                am_idx,exact_matches = ODBfilter.find_alias_matches(symbol,unfiltered_tsv,errors_fpath)
                am_df = unfiltered_tsv.loc[am_idx]
                em_df = ODBfilter.exact_match_df(unfiltered_tsv,exact_matches)
                if display_match_data:
                    print("alias_match")
                    display(am_df)
                    print("exact_match")
                    display(em_df)
                final_ksr_df = ODBfilter.select_known_species_records(symbol,em_df,am_df,ks_tids,unfiltered_fasta,
                                                                      manual_selections_fpath=tmp_manual_selections)
                if display_ksr:
                    print("Final known species records: ")
                    display(final_ksr_df)
                if symbol == 'CD151':
                    #Ensure non alias match/ exact match sequences not present in final ksr
                    self.assertFalse(len(final_ksr_df) == len(ks_tids))
                else:
                    self.assertTrue(len(final_ksr_df) == len(ks_tids))

                #Test manual selections cache
                if symbol == 'CALM1':
                    import contextlib
                    import io
                    print("Repeating CALM1 record selection. Should not ask for input (check this yourself)")

                    out_buf = io.StringIO()
                    print("Checking for cached selection output...")
                    with contextlib.redirect_stdout(out_buf):

                        final_ksr_df = ODBfilter.select_known_species_records(symbol, em_df, am_df, ks_tids,
                                                                              unfiltered_fasta,
                                                                              manual_selections_fpath=tmp_manual_selections)
                        cached_msg = 'To clear selections, either delete corresponding row in file at ' \
                                     '{0}'.format(tmp_manual_selections)
                        self.assertIn(cached_msg,out_buf.getvalue())
                    print("Cached selection output found.")

    def test_outgroup_selection(self):
        import SSfilter.ODBfilter
        test_symbol_list = ['ATP5MC1', 'CALM1', 'ATPIF1', 'CD151']
        tax_subset = ['10090_0', '43179_0', '9606_0', '10116_0', '42254_0', '9601_0']
        # symbol = 'ATP5MC1'
        symbol = 'IRF2BP2'
        errors_fpath = "{0}/outgroup_errors.tsv".format(test_tmp_dir)
        ks_tids = ['10090_0', '43179_0', '9606_0']
        tsv_inpath = "{0}/ODB/{1}.tsv".format(test_data_dir,symbol)
        unfiltered_tsv = SSfasta.load_tsv_table(tsv_inpath, tax_subset=tax_subset)
        unfiltered_fasta = "{0}/ODB/{1}.fasta".format(test_data_dir,symbol)
        am_idx, exact_matches = ODBfilter.find_alias_matches(symbol, unfiltered_tsv, errors_fpath)
        am_df = unfiltered_tsv.loc[am_idx]
        em_df = ODBfilter.exact_match_df(unfiltered_tsv, exact_matches)
        final_ksr_df = ODBfilter.select_known_species_records(symbol, em_df, am_df, ks_tids, unfiltered_fasta)

        final_dict = ODBfilter.select_outgrup_records(em_df,am_df,ks_tids,final_ksr_df,unfiltered_fasta)
        final_df = final_dict['final_df']
        assert(len(final_df) == len(tax_subset))

class NCBIFilterFunctionTest(unittest.TestCase):

    def test_NCBI_load(self):
        from SSfilter.NCBIfilter import load_NCBI_fasta_df

        test_fpath = "{0}/NCBI/9999/ATP5MC1.fasta".format(test_data_dir)
        taxid_dict = {'Urocitellus parryii':9999}
        ncbi_df = load_NCBI_fasta_df(test_fpath,taxid_dict)
        with pd.option_context('display.max_columns',None):
            self.assertTrue('XP_026242723.1' in ncbi_df.index)
            self.assertTrue(9999 in ncbi_df['organism_taxid'].unique())

    def test_select_NCBI_record(self):
        from SSfilter.NCBIfilter import select_NCBI_record
        from SSfilter.ODBfilter import process_ODB_input
        test_odb_path = "{0}/ODB/ATP5MC1.fasta".format(test_data_dir)
        test_ncbi_path = "{0}/NCBI/9999/ATP5MC1.fasta".format(test_data_dir)

        from SSutility import config, tax_subset, gene_id_df, tax_table
        tax_subset = ['10090_0', '43179_0', '9606_0', '10116_0', '42254_0', '9601_0']
        results = process_ODB_input('ATP5MC1',config,tax_subset)
        final_odb = results['final_df']
        taxid_dict = {config['NCBI']['NCBITaxName']:config['NCBI']['NCBITaxID']}
        final_combined = select_NCBI_record(test_odb_path,test_ncbi_path,taxid_dict,final_odb,['43179_0'])
        with pd.option_context('display.max_columns', None):
            # display(final_combined)
            self.assertTrue(len(final_combined) == len(final_odb)+1)
            self.assertTrue('Urocitellus parryii' in final_combined['organism_name'].unique())
            self.assertTrue('XP_026242723.1' in final_combined.index)

    def test_process_combined(self):
        from SSfilter.NCBIfilter import select_NCBI_record, combined_records_processing
        from SSfilter.ODBfilter import process_ODB_input
        from Bio import SeqIO,Seq

        test_symbol = "CALM1"

        test_odb_path = "{0}/ODB/{1}.fasta".format(test_data_dir,test_symbol)
        test_ncbi_path = "{0}/NCBI/9999/{1}.fasta".format(test_data_dir,test_symbol)

        config, tax_subset, gene_id_df, tax_table = SSconfig.config_initialization()
        tax_subset = ['10090_0', '43179_0', '9606_0', '10116_0', '42254_0', '9601_0']
        results = process_ODB_input(test_symbol, config, tax_subset)
        final_odb,em_df,am_df = results['final_df'],results['em_df'],results['am_df']
        taxid_dict = {config['NCBI']['NCBITaxName']: config['NCBI']['NCBITaxID']}

        final_combined = select_NCBI_record(test_odb_path, test_ncbi_path, taxid_dict, final_odb, ['43179_0'])

        test_unaln = "{0}/{1}_combined.fasta".format(test_tmp_dir,test_symbol)
        test_aln = "{0}/{1}_aligned.fasta".format(test_tmp_dir,test_symbol)
        test_records = "{0}/{1}_records.tsv".format(test_tmp_dir,test_symbol)

        comb_proc = combined_records_processing(config,am_df,em_df,final_combined,test_symbol,
                                    out_unaln_fasta=test_unaln,out_aln_fasta=test_aln,
                                    out_tsv_fpath=test_records)
        with pd.option_context('display.max_columns',None):
            display(comb_proc)
        #Check that final_combined record order is respected in outfile writing
        unaln_records = SeqIO.parse(test_unaln,'fasta')
        for i,fasta in enumerate(unaln_records):
            fc_idx_entry = final_combined.index[i]
            self.assertTrue(fc_idx_entry == fasta.id)
        aln_records = SeqIO.parse(test_aln, 'fasta')
        for i, fasta in enumerate(aln_records):
            fc_idx_entry = final_combined.index[i]
            self.assertTrue(fc_idx_entry == fasta.id)

        display_am = False
        if display_am:
            print("Alias matches")
            display(am_df.loc[:,['pub_gene_id','og_name','description']])


        self.assertTrue('manual selection' in comb_proc['selection_type'].unique())
        self.assertTrue(comb_proc.loc['43179_0:00184f','selection_type'] == "alias match min dist")
        self.assertTrue(comb_proc.loc['9601_0:003546', 'selection_type'] == "symbol match single record")


if __name__ == '__main__':

    from SSutility.SSdirectory import create_directory,empty_directory
    create_directory(test_tmp_dir)
    empty_directory(test_tmp_dir)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=[ImportWarning,DeprecationWarning])
        unittest.main(buffer=False)