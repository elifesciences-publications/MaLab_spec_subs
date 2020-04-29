#SSconfig.py - Funcitons for reading config files, including run parameters, gene symbols list, and species lists.
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

from SSutility import SSdirectory
import numpy as np
import pandas as pd
from IPython.display import display

def parse_config(config_file="config/config.txt"):
    """Parse config text file (INI format) to establish paramters for the run

    config_file: path to the config file ("config/config.txt" by default)
    """
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)
    return config


def parse_genes(genes_path="config/genes.txt"):
    """Parses gene file into list of uppercase, whitespace trimmed gene names.
    Currently deprecated - use gene_symbol column from read_geneID_file instead to get gene symbols.

    :param genes_path: file path to txt file list of gene symbols. Lines should be valid gene symbols or will
    be logged downstream as errors when attempting to fetch OrthoDB data.
    :return genes: list of gene symbols (capitalized and stripped of white space from genes_path)
    """
    gene_flines = open(genes_path).readlines()
    genes = [gene.strip().upper() for gene in gene_flines]
    return genes


def parse_species(species_path="config/v10_0_species.txt"):
    # Reads species list from file in config directory. Also returns a hashcode for the list of species used
    with open(species_path,'r') as spec_f:
        spec_lines = spec_f.readlines()
        spec_f.close()
    species = [spec.strip() for spec in spec_lines]
    concat = ""
    for spec in species:
        concat = concat + spec
    hc = np.abs(hash(concat))
    return species, hc

def read_geneID_file(csv_fpath):
    """Reads DataFrame from csv file stored at csv_path. The gene_symbol column will be used as the gene list for the
    run. human_gene_id will be used to fetch orthologs from the species specified by NCBITaxID in the config file
    (Urocitellus parryii by default but can be changed to fit user's needs). gene_symbol will be converted to upper case
    for file name consistency.

    :param csv_fpath: file path to table containing at least gene symbols and human gene IDs with mandatory column
    names gene_symbol and human_gene_id respectively
    :return: gene_id_df: DataFrame object containing fields from file
    """
    field_conv_dict = {"human_gene_id":str}
    gene_id_df = pd.read_csv(csv_fpath, dtype=field_conv_dict,index_col='overall_index')
    upper_symbol_srs = gene_id_df['gene_symbol'].str.upper()
    gene_id_df['gene_symbol'] = upper_symbol_srs
    return gene_id_df



def odb_tablev9(species_list, table_path="odb9v1_raw/odb9v1_species.tab"):
    """Reads orthodb v9 tsv file into a DataFrame of species names/ tax_ids and other ODB information.
        Mainly used for taxid <-> species name conversions
    """
    odb = pd.read_csv(table_path, delimiter="\t", header=None,
                      names=['tax_id", "odb_id", "spec_name", "clustered_genes", "ortho_groups", "mapping_type'])
    filtered = pd.DataFrame(columns=odb.columns)
    for spec in species_list:
        row = odb[odb['spec_name'] == spec]
        filtered = filtered.append(row)
    filtered.drop(columns=['clustered_genes", "ortho_groups", "mapping_type'], inplace=True)
    return filtered


def odb_tablev10(species_list, table_path="config/odb10v0_species.tab"):
    """odb10v0_species.tab:
    1.	NCBI tax id
    2.	Ortho DB individual organism id, based on NCBI tax id
    3.	scientific name inherited from the most relevant NCBI tax id
    4.	genome asssembly id, when available
    5.	total count of clustered genes in this species
    6.	total count of the OGs it participates
    7.	mapping type, clustered(C) or mapped(M)
    Reads above file into a DataFrame used for tax_id/ species name information, limited to species in species_list
    """
    odb = pd.read_csv(table_path, delimiter="\t", header=None,
                      names=['tax_id', 'odb_id', 'spec_name', 'assembly_id', 'clustered_genes', 'ortho_groups',
                             'mapping_type'])
    filtered = pd.DataFrame(columns=odb.columns)
    for spec in species_list:
        row = odb[odb['spec_name'] == spec]
        filtered = filtered.append(row)
    filtered.drop(columns=['clustered_genes', 'ortho_groups', 'mapping_type'], inplace=True)
    filtered = filtered.sort_values(by='tax_id')
    return filtered

def config_initialization():
# Read config files
    config = parse_config()
    run_config = config['RUN']
    tax_subset = list(config['AnalysisODBTaxSubset'].keys())
    gene_id_fpath = run_config['IDFilePath']
    species_path = run_config['SpeciesFilePath']
    spec_list, hc = parse_species(species_path)
    gene_id_df = read_geneID_file(gene_id_fpath)
    tax_table = odb_tablev10(spec_list)
    run_name = run_config['RunName']
    SSdirectory.create_run_directory(run_name)
    #Copy config file to summary_directory
    src_config = "config/config.txt"
    dest_config = "{0}/summary/config.txt".format(run_name)
    import shutil
    shutil.copy(src_config,dest_config)
    return config, tax_subset, gene_id_df, tax_table

def main():
    DISPLAY_PARAMS = False
    if DISPLAY_PARAMS:
        config, spec_list, tax_subset, gene_id_df, tax_table = config_initialization()
        run_config, odb_config = config['RUN'], config['ODB']
        run_name = run_config['RunName']
        test_species = odb_config['ODBTestSpecies']
        species_path = run_config['SpeciesFilePath']
        print("Tax table for species list at {0}".format(species_path))
        with pd.option_context("display.max_columns", None):
            display(tax_table)
        print("Analysis taxonomy subset: {0}".format(list(tax_subset.keys())))
        print("Gene ID table")
        # print(gene_id_df['gene_symbol'].values)
        display(gene_id_df)
        print("Run Name: " + run_name)
        # Verify that species table, gene list, and run_name are correct

if __name__ == '__main__':
    main()