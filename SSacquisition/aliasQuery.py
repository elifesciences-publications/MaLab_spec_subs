import numpy as np
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import re
import os
from SSutility.SSerrors import GeneCardsError, write_errors, load_errors, print_errors

def write_aliases_f(aliases, aliases_fpath):
    """Write aliases data from GeneCards to a txt file list at aliases_fpath"""
    aliases_f = open(aliases_fpath, 'wt')
    for a in aliases:
        aliases_f.write(a.strip() + '\n')
    aliases_f.close()


def alias_GC_query(driver,gene_name):
    """Queries GeneCards for alias data for gene_name. Should only be called if aliases_fpath doesn't exist
    (ie if query has not been previously run and written to file). Attempts GeneCards query - if gene_name
    leads to a single entry page, pulls aliases from page html. If query leads to a query results page,
    checks all linked entries to see if any contain gene_name. If none do (or other WebDriver issues arise),
    raises a GeneCardsError.

    If query was successful (either single result page or successfully chose linked result from query results),
    return aliases and gc_name (the gene identifier used by GeneCards). gc_name stored separately since alias html
    extraction will miss it otherwise. Also writes alias data to aliases_fpath

    Updated 01/10/2020. If function is consistently failing, check xpath class names against orthodb website

    :param driver: Selenium webdriver object; intialized previously with headless option
    :param gene_name: gene symbol which will be queried via GeneCards website for alias names
    :return: aliases: list of GeneCards alias names for gene_name
    :return gc_name: primary alias for this symbol used by GeneCards """


    gene_cards_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={0}".format(gene_name.upper())
    aliases_fpath = "alias_data/{0}_aliases.txt".format(gene_name)
    list_xpath = "//ul[@class='list-unstyled list-spacious']/li"
    elem_xpaths = [list_xpath]
    driver.get(gene_cards_url)
    aliases = []
    for xpath in elem_xpaths:
        elems = driver.find_elements_by_xpath(xpath)
        innerHTMLs = [elem.get_attribute("innerHTML") for elem in elems]
        col_aliases = [BeautifulSoup(markup,features="html5lib").find(text=True).strip() for markup in innerHTMLs]
        aliases.extend(col_aliases)
    if len(aliases) > 0:
        # Means gene_name query to GeneCards autoredirected to a single page - normal aliases scraping
        # HTML parsing for GeneCards website - end result is list of trimmed alias strings
        gc_re = re.search("gene=([A-Z0-9\-]+)", driver.current_url)
        gc_name = gc_re.groups()[0].strip()
        if gc_name not in aliases:
            aliases.insert(0, gc_name)
        # Cache aliases to aliases_fpath
        write_aliases_f(aliases, aliases_fpath)
    else:
        # Try search results page for gene_name; raise GeneCardsError if no results or check each page
        # for alias matching gene_name otherwise
        query_url = "https://www.genecards.org/Search/Keyword?queryString={0}".format(gene_name)
        driver.get(query_url)
        links_xpath = "//td[@class='gc-gene-symbol gc-highlight symbol-col']/a"
        link_elems = driver.find_elements_by_xpath(links_xpath)
        link_hrefs = [elem.get_attribute("href") for elem in link_elems]
        if link_elems:
            for elem_href in link_hrefs:
                driver.get(elem_href)
                link_url = driver.current_url
                elem_gc_name = re.search("gene=([A-Z0-9\-]+)", link_url).groups()[0].strip()
                elem_aliases = []
                for xpath in elem_xpaths:
                    elems = driver.find_elements_by_xpath(xpath)
                    innerHTMLs = [elem.get_attribute("innerHTML") for elem in elems]
                    col_aliases = [BeautifulSoup(markup,features="html5lib").find(text=True).strip() for markup in innerHTMLs]
                    elem_aliases.extend(col_aliases)
                if gene_name in elem_aliases or gene_name == elem_gc_name:
                    # Found query result with gene_name
                    if elem_gc_name not in elem_aliases:
                        elem_aliases.insert(0, elem_gc_name)
                    write_aliases_f(elem_aliases, aliases_fpath)
                    return
        # If either no link_elems (empty search results page), or none correspond to gene_name:
        raise GeneCardsError(0, "Could not automatically fetch alias data from GeneCards - consider searching manually")

def download_alias_data(gene_list, config):
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    window_size = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % window_size)
    driver = webdriver.Chrome(chrome_options=chrome_options)
    # driver.implicitly_wait(5)
    run_name,errors_fname = config['RUN']['RunName'],config['RUN']['ErrorsFileName']
    errors_fpath = "{0}/{1}".format(run_name,errors_fname)
    check_error_file, gc_errors_df = load_errors(errors_fpath,"GeneCardsError")
    for gene_name in gene_list:
        aliases_fpath = "alias_data/{0}_aliases.txt".format(gene_name)
        if not os.path.exists(aliases_fpath):
            if check_error_file and gene_name in gc_errors_df["gene_symbol"].unique():
                print_errors(gc_errors_df, gene_name)
            else:
                try:
                    alias_GC_query(driver, gene_name)
                except GeneCardsError as gc_error:
                    write_errors(errors_fpath,gene_name,gc_error)
    driver.quit()

def __clean_alias_data_dir(src_dir_path, target_dir_path, gene_symbols):
    """Creates a new directory at target_dir_path, copies any _alias.txt files from src_dir_path to target_dir_path
    if the corresponding gene symbol is in gene_symbols. Will not create directory if target_dir_path exists.

    :param target_dir_path:
    :param gene_symbols: iterable, must support check if gene symbol in gene_symbols (must provide .unique() for Pandas
    Series)
    :return count: number of files from src_dir_path copied to target_dir_path
    """
    try:
        os.mkdir(target_dir_path)
    except FileExistsError as fe:
        print("target_dir_path {0} already exists. Provide different target directory".format(target_dir_path))
        return
    count = 0
    for f_name in os.listdir(src_dir_path):
        f_symbol_match = re.search("(\w+)_aliases.txt", f_name)
        if f_symbol_match:
            symbol = f_symbol_match.groups()[0]
            if symbol in gene_symbols:
                count += 1
                src_fpath = os.path.join(src_dir_path, f_name)
                target_fpath = os.path.join(target_dir_path, f_name)
                if not os.path.exists(target_fpath):
                    os.rename(src_fpath, target_fpath)
    return count

