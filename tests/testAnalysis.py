import os
os.chdir("..")
import unittest
from SSutility import SSfasta
from IPython.display import display

class testCalculations(unittest.TestCase):

    def test_unique_pos(self):
        from SSanalysis.SSanalysiscalc import find_uniques
        test_outdir = "tests/test_data/output/ATP5MC1/"
        test_records = "{0}/ATP5MC1_records.tsv".format(test_outdir)
        test_msa = "{0}/ATP5MC1_msa.fasta".format(test_outdir)
        align_df = SSfasta.align_fasta_to_df(test_msa)

        filtered = align_df.iloc[:-1,:]

        sub_freq = 1
        test_idx = '43179_0:00103c'
        uniques = find_uniques(filtered,sub_freq,test_idx,False)
        for expr in [16 in uniques.columns,uniques.loc[test_idx,33]=='L',4 not in uniques.columns]:
            self.assertTrue(expr)

    # @unittest.skip("Dependency on Capra and Singh 2007 SSanalysis/score_conservation.py")
    def test_sc_calculations(self):
        from SSanalysis import score_conservation
        from SSanalysis import aas, blosum62_bg, sim_matrix

        test_outdir = "tests/test_data/output/ATP5MC1/"
        test_msa = "{0}/ATP5MC1_msa.fasta".format(test_outdir)
        align_df = SSfasta.align_fasta_to_df(test_msa)
        msa_nd = align_df.values
        sc_sw = score_conservation.calculate_sequence_weights(msa_nd)

        from SSanalysis.JSDcalc import calculate_sequence_weights, JSD, weighted_fcpseudo
        native_sw = calculate_sequence_weights(msa_nd)

        for i in range(len(native_sw)):
            self.assertAlmostEqual(native_sw[i],sc_sw[i])


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


    def test_SScalcs(self):
        from SSanalysis.SSanalysiscalc import generate_jsd_series

        test_outdir = "tests/test_data/output/ATP5MC1/"
        test_msa = "{0}/ATP5MC1_msa.fasta".format(test_outdir)
        align_df = SSfasta.align_fasta_to_df(test_msa)
        align_df = align_df.drop(index="XP_026242723.1")
        msa_nd = align_df.values
        test_idx = "43179_0:00103c"
        jsd,jsd_z = generate_jsd_series(test_idx,align_df)
        display(jsd.loc[30:50])
        display(align_df.loc[:,30:50])



if __name__ == "__main__":
    unittest.main()