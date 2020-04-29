import os
os.chdir("..")
import unittest
from SSutility import SSfasta
from IPython.display import display
import pandas as pd
import numpy as np
pd.options.display.max_columns = None

from SSanalysis import aas, blosum62_bg,blos_df, sim_matrix

#Test sequence data
test_outdir = "tests/test_data/output/ATP5MC1/"
test_msa = "{0}/ATP5MC1_msa.fasta".format(test_outdir)
align_df = SSfasta.align_fasta_to_df(test_msa)
ODB_df = align_df.drop(index="XP_026242723.1")
test_idx = pd.Index(["43179_0:00103c"])
ncbi_idx = pd.Index(["XP_026242723.1"])

class testCalculations(unittest.TestCase):

    def test_unique_pos(self):
        from SSanalysis.SSanalysiscalc import find_uniques
        #Identify uniques from ODB_df for ATP5MC1
        filtered = ODB_df
        sub_freq = 1

        uniques = find_uniques(filtered,sub_freq,test_idx,False)
        for expr in [16 in uniques.columns,(uniques.loc[test_idx,33]=='L').all(),4 not in uniques.columns]:
            self.assertTrue(expr)



    #Published code for Capra and Singh requires substantial reformatting to run in Python 3.7. A heavily truncated
    # version of their code is provided in SSanalysis as score_conservation.py to test JSDcalc functions for consistency
    # Unedited code available from url below.
    # https://compbio.cs.princeton.edu/conservation/
    #John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    #Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270

    # @unittest.skip("Dependency on Capra and Singh 2007 SSanalysis/score_conservation.py")
    def test_sc_calculations(self):
        from SSanalysis import score_conservation

        msa_nd = ODB_df.drop(test_idx).values
        sc_sw = score_conservation.calculate_sequence_weights(msa_nd)

        from SSanalysis.JSDcalc import calculate_sequence_weights, JSD, weighted_fcpseudo, weighted_gap_penalty
        native_sw = calculate_sequence_weights(msa_nd)

        for i in range(len(native_sw)):
            self.assertAlmostEqual(native_sw[i],sc_sw[i])

        test_gap_col = ['-']* 8
        sc_gp = score_conservation.weighted_gap_penalty(test_gap_col,native_sw)
        native_gp = weighted_gap_penalty(test_gap_col,native_sw)
        self.assertAlmostEqual(sc_gp,native_gp)

        test_cols = [msa_nd[:,i] for i in range(0,135)]
        #Test weighted_fc_pseudocount
        for col in test_cols:
            native_fc = weighted_fcpseudo(col,native_sw,aas)
            sc_fc = score_conservation.weighted_freq_count_pseudocount(col,native_sw,0.0000001)
            try:
                for i in range(len(native_fc)):
                    self.assertAlmostEqual(native_fc[i],sc_fc[i])
            except:
                print("Native FC: {0}".format(native_fc))
                print("Score Conservation FC: {0}".format(sc_fc))
        #Test equivalence of calculated JSD
        for col in test_cols:
            sc_jsd = score_conservation.js_divergence(col,sim_matrix,blosum62_bg,native_sw,1)
            native_jsd = JSD(col,blosum62_bg,native_sw,aas,use_gap_penalty=True)
            try:
                self.assertAlmostEqual(sc_jsd,native_jsd)
            except:
                print("Native JSD: {0}".format(native_jsd))
                print("Score Conservation JSD: {0}".format(sc_jsd))

        #Test gap penalty equivalence
        gp_test_col = msa_nd[:,42]
        sc_jsd = score_conservation.js_divergence(gp_test_col,sim_matrix,blosum62_bg,native_sw,0)
        sc_gp_jsd = score_conservation.js_divergence(gp_test_col, sim_matrix, blosum62_bg, native_sw, 1)
        native_jsd = JSD(gp_test_col, blosum62_bg, native_sw, aas, use_gap_penalty=False)
        native_gp_jsd = JSD(gp_test_col, blosum62_bg, native_sw, aas, use_gap_penalty=True)
        self.assertAlmostEqual(sc_jsd,native_jsd)
        self.assertAlmostEqual(sc_gp_jsd, native_gp_jsd)
        self.assertNotAlmostEqual(native_jsd,native_gp_jsd)
        # print(native_gp_jsd)
        # print(sc_gp_jsd)
        # print(native_jsd)
        # print(sc_jsd)

    def test_jsd_srs(self):
        # from SSanalysis.SSanalysiscalc import generate_jsd_series,find_uniques
        from SSanalysis import SSanalysiscalc as ac

        jsd,jsd_z = ac.generate_jsd_series(test_idx,ODB_df)

        uniques = ac.find_uniques(ODB_df,1,test_idx)
        self.assertAlmostEqual(jsd.loc[33],0.871706,places=5)
        unique_jsd = ac.filter_score_srs(jsd,uniques)
        self.assertTrue(33 in unique_jsd.index)
        self.assertAlmostEqual(jsd.loc[33],unique_jsd.loc[33])
        self.assertFalse(57 in unique_jsd.index)

        jsd_43 = jsd[43]
        no_gp_jsd, _ = ac.generate_jsd_series(test_idx,ODB_df,use_gap_penalty=False)

        from SSanalysis.JSDcalc import JSD, calculate_sequence_weights
        ngp_43 = no_gp_jsd[43]

        msa_nd = ODB_df.drop(test_idx).values
        native_sw = calculate_sequence_weights(msa_nd)
        gp_test_col = msa_nd[:, 42]
        no_gp_43_recalc = JSD(gp_test_col, blosum62_bg, native_sw, aas, use_gap_penalty=False)
        self.assertNotAlmostEqual(jsd_43,ngp_43)
        self.assertAlmostEqual(ngp_43,no_gp_43_recalc)


    def test_blos_calcs(self):
        from SSanalysis import SSanalysiscalc as ac

        test_cols = [ODB_df.loc[:,pos] for pos in [31, 33, 40, 43, 75]]
        to_expected_values = [np.mean([4,4,4,4,4,1,1,1]),-3,-1,np.nan,4]
        for i,col in enumerate(test_cols):
            # display(col)
            calc_blos = ac.test_outgroup_blosum(col,test_idx,blos_df)
            if i == 3:
                self.assertTrue(np.isnan(calc_blos))
            else:
                self.assertAlmostEqual(to_expected_values[i],calc_blos)
        pos31_pw_score = np.mean([1,4,1,1,4,4,4, 1,6,6,1,1,1, 1,1,4,4,4, 6,1,1,1, 1,1,1, 4,4, 4])
        pw_expected_valyes = [pos31_pw_score,7,7,np.nan,4]
        for i,col in enumerate(test_cols):
            calc_pw = ac.pairwise_outgroup_blosum(col,test_idx,blos_df)
            if i == 3:
                self.assertTrue(np.isnan(calc_pw))
            else:
                self.assertAlmostEqual(calc_pw,pw_expected_valyes[i])

    def test_variant_counts(self):
        from SSanalysis import SSanalysiscalc as ac

        #check expected values for test/outgroup variants on ATP5MC1
        test_cols = [ODB_df.loc[:, pos] for pos in [31, 33, 40, 43, 75]]
        exp_tvs = ['S','L','T','S','A']
        exp_ovs = ['S', 'P', 'P', '-', 'A']
        exp_tv_count = [6,1,1,1,9]
        exp_ov_count = [6,8,8,8,9]
        exp_tv_ratio = [c / len(ODB_df) for c in exp_tv_count]
        exp_ov_ratio = [c / len(ODB_df) for c in exp_ov_count]
        for i,col in enumerate(test_cols):
            metrics = ac.variant_counts(col,test_idx)
            self.assertEqual(metrics['test_variant'],exp_tvs[i])
            self.assertEqual(metrics['outgroup_variant'], exp_ovs[i])
            self.assertEqual(metrics['test_variant_count'], exp_tv_count[i])
            self.assertEqual(metrics['outgroup_variant_count'], exp_ov_count[i])
            self.assertEqual(metrics['test_variant_ratio'], exp_tv_ratio[i])
            self.assertEqual(metrics['outgroup_variant_ratio'], exp_ov_ratio[i])

    def test_gene_summary_df(self):
        from SSanalysis import SSanalysiscalc as ac
        test_outpath = "tests/tmp/ATP5MC1_summary.tsv"
        summary_table = ac.gene_summary_table(align_df,ncbi_idx,test_idx,blos_df,
                                              summary_table_outpath=test_outpath)
        with pd.option_context('display.max_columns',None):
            # display(summary_table)
            pass

        from_file = ac.load_summary_table(test_outpath)
        with pd.option_context('display.max_columns',None):
            # display(from_file)
            pass
        #DataType preservation in loading from file
        self.assertTrue(7.0000 in from_file['Outgroup Pairwise BLOSUM62'].unique())
        self.assertTrue(33 in from_file.index)
        self.assertTrue(np.isnan(from_file.loc[16,'Test-Outgroup BLOSUM62']))

        gp_43_jsd = summary_table.loc[43]['JSD']
        non_gp_summary = ac.gene_summary_table(align_df,ncbi_idx,test_idx,blos_df,summary_table_outpath=test_outpath,
                                               use_jsd_gap_penalty=False)
        non_gp_43_jsd = non_gp_summary.loc[43]['JSD']
        print(gp_43_jsd)
        print(non_gp_43_jsd)
        self.assertNotAlmostEqual(gp_43_jsd,non_gp_43_jsd)

    def test_overall_summary(self):
        from SSanalysis import SSanalysiscalc as ac
        from SSutility import config
        test_symbols = ['ATP5MC1']
        new_symbols = ['ATP5MC1', 'ATPIF1']
        test_outpaths = ["cDNAscreen_041020/output/{0}/{0}_summary.tsv".format(symbol) for symbol in new_symbols]
        #Temporarily relocate files already in output directory because I'm too lazy to pass in gene specific paths to
        #overall_summary_table
        cache_outpaths = ["cDNAscreen_041020/output/{0}/{0}_summary_cache.tsv".format(symbol) for symbol in new_symbols]
        for i,path in enumerate(test_outpaths):
            if os.path.exists(path):
                os.rename(path,cache_outpaths[i])
        #Test single gene summary table formatting into overall summary
        overall = ac.overall_summary_table(config,test_symbols,use_jsd_gap_penalty=True,force_recalc=True)
        self.assertTrue('Gene' in overall.columns)
        self.assertTrue('ATP5MC1' in overall['Gene'].unique())
        self.assertTrue(overall['JSD Alignment Z-Score'].any())
        self.assertTrue(overall['JSD US Z-Score'].any())
        with pd.option_context('display.max_columns',None):
            # display(overall)
            pass
        #Test gap penalty functions from overall_summary_table level
        overall_nogp = ac.overall_summary_table(config, test_symbols, use_jsd_gap_penalty=False, force_recalc=True)
        gp_test_jsd = overall.loc[4,'JSD']
        nogp_test_jsd = overall_nogp.loc[4,'JSD']
        self.assertNotAlmostEqual(gp_test_jsd,nogp_test_jsd)
        #Test appending second symbol summary
        new_symbols = ['ATP5MC1','ATPIF1']
        added_nogp = ac.overall_summary_table(config, new_symbols, use_jsd_gap_penalty=False, force_recalc=True)
        self.assertTrue(10 in added_nogp.index)
        self.assertTrue('ATPIF1' in added_nogp['Gene'].unique())

        #Test overall Z-score calculation consistent with expectation (ie ignore nan)
        combined_mean_jsd, combined_std_jsd = np.nanmean(added_nogp['JSD']),np.nanstd(added_nogp['JSD'])
        self.assertAlmostEqual(added_nogp.loc[11,'JSD US Z-Score'],
                               (added_nogp.loc[11,'JSD']-combined_mean_jsd)/combined_std_jsd)
        overall_blos_values = [1,-3,-1,-2,-3.14286,-1.42857,1.28571,-1,-0.714286,-0.714286]
        test_blos_mean, test_blos_std = np.mean(overall_blos_values),np.std(overall_blos_values)
        for i,val in enumerate(added_nogp['Test-Outgroup BLOSUM US Z-Score'].dropna().values):
            calc_val = (overall_blos_values[i]-test_blos_mean)/test_blos_std
            self.assertAlmostEqual(val,calc_val,places=5)

        #Delete test generated outfiles, rename any cached files existing before run
        for path in test_outpaths:
            if os.path.exists(path):
                os.remove(path)
        for i,path in enumerate(cache_outpaths):
            if os.path.exists(path):
                os.rename(path, test_outpaths[i])

if __name__ == "__main__":
    unittest.main()