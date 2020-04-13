import os
import glob
import shutil

"""File and directory management functions"""


def create_directory(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def empty_directory(path):
    for i in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(i):
            shutil.rmtree(i)
        else:
            os.remove(i)


def create_run_directory(run_name):
    """Make diretory tree for a run. Also creates tmp directory for storage of temporary alignment files"""
    dirpaths = ["{0}", "{0}/input", "{0}/input/ODB", "{0}/input/NCBI", "{0}/output", "{0}/summary", "{0}/run_params"]
    for dirpath in dirpaths:
        formatted = dirpath.format(run_name)
        create_directory(formatted)
    create_directory("tmp")
    empty_directory("tmp")


def write_run_params_file(config, spec_path, spec_hc):
    """Documents some run specific parameters.

    config: config file object storing some run specifications (directory names, file paths)
    spec_path: file path for the input species list being used for analysis
    spec_hc: hashcode generated from species list.
    """
    config_keys = ["RunName", "GenesFilePath", "ODBLevel"]
    run_name = config["RunName"]
    fpath = run_name + "/run_params/params.txt"
    params_f = open(fpath, 'wt')
    for key in config_keys:
        val = config[key]
        file_line = "{0}: {1}\n".format(key, val)
        params_f.write(file_line)
    params_f.write("species_list: {0}\n".format(spec_path))
    params_f.write("species_hashcode: {0}\n".format(spec_hc))
    params_f.close()

def convert_fname_uppercase(file_dir,pat):
    """Converts file names in file_dir to upper case gene symbols (file extension will remain intact).

    :param file_dir: directory in which file names will be converted to upper case
    :param pat: regular expression string which has two groups, the gene symbol and the rest of the file name
    :return: None
    """
    import re
    for fname in os.listdir(file_dir):
        match = re.search(pat,fname)
        if match:
            symbol = match.groups()[0]
            file_ext = match.groups()[1]
            if not symbol.upper() == symbol:
                src_fpath = os.path.join(file_dir,fname)
                target_fpath = os.path.join(file_dir,symbol.upper()+file_ext)
                os.replace(src_fpath,target_fpath)


# importlib.reload(SSdirectory)
#run_name = "cDNAscreen_041020"
# NCBI_pat = "(\w+)(\.fasta)"
# ODB_pat = "(\w+)(\.fasta|\.tsv)"
# SSdirectory.convert_fname_uppercase("{0}/input/NCBI/9999".format(run_name),NCBI_pat)
# SSdirectory.convert_fname_uppercase("{0}/input/ODB".format(run_name),ODB_pat)