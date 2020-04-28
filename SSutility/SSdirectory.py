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
    create_directory("alias_data")


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


def main():
    run_name = "cDNAscreen_041020"
    NCBI_pat = r"(\w+)(\.fasta)"
    ODB_pat = r"(\w+)(\.fasta|\.tsv)"
    convert_fname_uppercase("{0}/input/NCBI/9999".format(run_name),NCBI_pat)
    convert_fname_uppercase("{0}/input/ODB".format(run_name),ODB_pat)

if __name__ == '__main__':
    main()