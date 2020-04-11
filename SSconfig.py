# Funcitons for reading config files, including run parameters, gene symbols list, and species lists.
import SSdirectory
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
    return config["DEFAULT"]


def parse_genes(genes_path="config/genes.txt"):
    """Parses gene file into list of uppercase, whitespace trimmed gene names"""
    gene_flines = open(genes_path).readlines()
    genes = [gene.strip().upper() for gene in gene_flines]
    return genes


def parse_species(species_path="config/v10_0_species.txt"):
    # Reads species list from file in config directory. Also returns a hashcode for the list of species used
    spec_lines = open(species_path).readlines()
    species = [spec.strip() for spec in spec_lines]
    concat = ""
    for spec in species:
        concat = concat + spec
    hc = np.abs(hash(concat))
    return species, hc


def odb_tablev9(species_list, table_path="odb9v1_raw/odb9v1_species.tab"):
    """Reads orthodb v9 tsv file into a DataFrame of species names/ tax_ids and other ODB information.
        Mainly used for taxid <-> species name conversions
    """
    odb = pd.read_csv(table_path, delimiter="\t", header=None,
                      names=["tax_id", "odb_id", "spec_name", "clustered_genes", "ortho_groups", "mapping_type"])
    filtered = pd.DataFrame(columns=odb.columns)
    for spec in species_list:
        row = odb[odb["spec_name"] == spec]
        filtered = filtered.append(row)
    filtered.drop(columns=["clustered_genes", "ortho_groups", "mapping_type"], inplace=True)
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
    Reads above file into a DataFrame used for tax_id/ species name information
    """
    odb = pd.read_csv(table_path, delimiter="\t", header=None,
                      names=["tax_id", "odb_id", "spec_name", "assembly_id", "clustered_genes", "ortho_groups",
                             "mapping_type"])
    filtered = pd.DataFrame(columns=odb.columns)
    for spec in species_list:
        row = odb[odb["spec_name"] == spec]
        filtered = filtered.append(row)
    filtered.drop(columns=["clustered_genes", "ortho_groups", "mapping_type"], inplace=True)
    return filtered

def config_initialization():
# Read config files
    config = parse_config()

    species_path = "config/v10_0_species.txt"
    spec_list, hc = parse_species(species_path)
    gene_list = parse_genes(config["GenesFilePath"])
    tax_table = odb_tablev10(spec_list)
    run_name = config["RunName"]
    SSdirectory.create_run_directory(run_name)
    return config, spec_list, gene_list, tax_table

DISPLAY_PARAMS = False
if DISPLAY_PARAMS:
    config, spec_list, gene_list, tax_table = config_initialization()
    run_name = config["RunName"]
    test_species = config["ODBTestSpecies"]
    species_path = config["SpeciesFilePath"]
    print("Tax table for species list at {0}".format(species_path))
    with pd.option_context("display.max_columns", None):
        display(tax_table)
    print("Gene list: " + str(gene_list))
    print("Run Name: " + run_name)
    # Verify that species table, gene list, and run_name are correct