# spec_subs_main.py - Handles data acquisition, filtering, substitution analysis, and plotting for
#comparison of 13LGS and AGS specific protein sequence changes.
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

from SSutility import SSdirectory, SSerrors
import sys, os

def ss_acquisition(config, gene_id_df, tax_table):
    """Uses gene_symbol column in gene_id_df to query OrthoDB and NCBI for input sequence data. Downloads alias
    information for all symbols from GeneCards into alias_data.

    :param config: configparser object with run params. See SSutility.SSconfig
    :param gene_id_df: DataFrame containing Gene ID and symbol information. See SSutility.SSconfig
    :param tax_table: OrthoDB taxonomy table filtered down to input species list in config[SpeciesFilePath].
    See SSutility.SSconfig
    :return: N/A
    """
    from SSacquisition import ODBquery,NCBIquery,aliasQuery
    print("===Record Data Acquisition====")
    gene_symbols = gene_id_df['gene_symbol']
    valid_queries, failed_queries = ODBquery.download_ODB_input(gene_symbols, tax_table, config)
    ags_mapped_df = NCBIquery.download_AGS_data(gene_id_df, config)
    aliasQuery.download_alias_data(gene_symbols, config)

def ss_filter(config, tax_subset, gene_id_df):
    """Record Filtering for all input data. Logs error and quality check output, generates final dataset fastas and
    records.tsv in [run_dir]/output, skipping OrthoDB/NCBI query errors and files with no usable test species data.

    :param config: configparser object, see SSutility.SSconfig
    :param tax_subset: list of Taxonomy IDs to be included in filtered output, read from config.txt
    [AnalysisODBTaxSubset] section. Can be configured to select amy subset of species from raw OrthoDB input.
    :param gene_id_df: DataFrame with gene ID and symbol information, SSutility.SSconfig
    :return: N/A. Writes filtered output files to [run_dir]/output subdirectories separated by gene symbol.
    """

    from SSfilter import ODBfilter,NCBIfilter
    print("===Record Data Filtering===")
    print("Sequence Data Error and Quality Check warnings for mouse/human/13LGS will be printed")
    print("Error logs and QC can be found within run directory ")

    gene_symbols = gene_id_df['gene_symbol']
    run_name,errors_fname,qc_fname= [config['RUN'][key] for key in ['RunName','ErrorsFileName','QCFileName']]
    errors_fpath = "{0}/{1}".format(run_name,errors_fname)
    qc_fpath = "{0}/{1}".format(run_name,qc_fname)

    check_errors,errors_df = SSerrors.load_errors(errors_fpath)
    check_qc, qc_df = SSerrors.load_errors(qc_fpath)
    qc_symbols = []

    for symbol in gene_symbols:
        out_records_fpath = "{0}/output/{1}/{1}_records.tsv".format(run_name, symbol)
        if not os.path.exists(out_records_fpath) or config['RUN'].getboolean('OverwriteFilter'):
            if check_errors and symbol in errors_df['gene_symbol'].unique():
                odb_error = (symbol in errors_df.loc[errors_df['error_type']=="OrthoDBQueryError",'gene_symbol'].unique())
                ncbi_error = (symbol in errors_df.loc[errors_df['error_type'] == "NCBIQueryError", 'gene_symbol'].unique())
                if odb_error or ncbi_error:
                    continue
                elif symbol in errors_df.loc[errors_df['error_type']=="SequenceDataError",'gene_symbol'].unique():
                    SSerrors.print_errors(errors_df,symbol,error_type="SequenceDataError")
                    continue
            try:
                NCBIfilter.final_combined_input(config,symbol,tax_subset)
            except SSerrors.SequenceDataError as e:
                continue
        #Write QC output for previously filtered results (ie out_records_fpath exists)
        elif check_qc and symbol in qc_df['gene_symbol'].unique():
            qc_symbols.append(symbol)
    print("==Previosuly cached QC==")
    for qc_symbol in qc_symbols:
        for i,row in qc_df.loc[qc_df['gene_symbol']==qc_symbol,:].iterrows():
            print("{0}\t{1}".format(row.values[0],row.values[1]))

def ss_analysis(config,gene_id_df):
    """Calculates jensen-shannon divergence, BLOSUM62 scores, and various other alignment metrics for the filtered
    dataset. JSD calculation code is modified from the below paper and accompanying code and is found in
    SSanalysis/JSDcalc.py

    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :return:
    """
    from SSanalysis import SSanalysiscalc as ac
    gene_symbols = gene_id_df['gene_symbol']
    ac.overall_summary_table(config, gene_symbols, use_jsd_gap_penalty=True, force_recalc=False)

def ss_visualization(config):
    from SSanalysis import SSvisualization
    plot_symbols = ["ATP5MC1", "ATP5MC3", "ATPIF1", "MANF"]
    SSvisualization.summary_plots(config,plot_symbols)


def main():
    import SSutility
    from SSutility import config, tax_subset, gene_id_df, tax_table
    ss_acquisition(config, gene_id_df, tax_table)
    ss_filter(config, tax_subset, gene_id_df)
    ss_analysis(config,gene_id_df)
    ss_visualization(config)

if __name__ == '__main__':
    main()
