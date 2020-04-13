import numpy as np
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import re
from SSerrors import GeneCardsError, write_errors

def format_odb_field(field):
    """Remove spaces, commas, and capitalization from alias/ odb fields to search for string matches.
    If field is empty (np.nan), return empty string"""
    if (type(field)) == str:
        field = field.replace(" ", "")
        field = field.replace(",", "")
        field = field.replace("\n", "")
        return field.lower()
    elif type(field) == float and np.isnan(field):
        return ""


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
    aliases_fpath = "aliases_data/" + gene_name + "_aliases.txt"

    list_xpath = "//ul[@class='list-unstyled list-spacious']/li"
    elem_xpaths = [list_xpath]
    driver.get(gene_cards_url)
    aliases = []
    for xpath in elem_xpaths:
        elems = driver.find_elements_by_xpath(xpath)
        innerHTMLs = [elem.get_attribute("innerHTML") for elem in elems]
        col_aliases = [BeautifulSoup(markup).find(text=True).strip() for markup in innerHTMLs]
        aliases.extend(col_aliases)
    if len(aliases) > 0:
        # Means gene_name query to GeneCards autoredirected to a single page - normal aliases scraping
        # HTML parsing for GeneCards website - end result is list of trimmed alias strings
        gc_re = re.search("gene=([A-Z0-9]+)", gene_cards_url)
        gc_name = gc_re.groups()[0].strip()
        if gc_name not in aliases:
            aliases.insert(0, gc_name)
        # Cache aliases to aliases_fpath
        write_aliases_f(aliases, aliases_fpath)
        driver.quit()
    else:
        # Try search results page for gene_name; raise GeneCardsError if no results or check each page
        # for alias matching gene_name otherwise
        query_url = "https://www.genecards.org/Search/Keyword?queryString={0}".format(gene_name)
        driver.get(query_url)
        links_xpath = "//td[@class='gc-gene-symbol gc-highlight symbol-col']/a"
        link_elems = driver.find_elements_by_xpath(links_xpath)
        if link_elems:
            for elem in link_elems:
                elem_href = elem.get_attribute("href")
                driver.get(elem_href)
                query_url = driver.current_url
                elem_gc_name = re.search("gene=([A-Z0-9]+)", query_url).groups()[0].strip()
                elem_aliases = []
                for xpath in elem_xpaths:
                    elems = driver.find_elements_by_xpath(xpath)
                    innerHTMLs = [elem.get_attribute("innerHTML") for elem in elems]
                    col_aliases = [BeautifulSoup(markup).find(text=True).strip() for markup in innerHTMLs]
                    elem_aliases.extend(col_aliases)
                if gene_name in elem_aliases or gene_name == elem_gc_name:
                    # Found query result with gene_name
                    driver.quit()
                    if elem_gc_name not in elem_aliases:
                        elem_aliases.insert(0, elem_gc_name)
                    write_aliases_f(elem_aliases, aliases_fpath)
                    return
        # If either no link_elems (empty search results page), or none correspond to gene_name:
        driver.quit()
        raise GeneCardsError(0, "Could not automatically fetch alias data from GeneCards - consider searching manually")

def download_alias_data(gene_list, config):
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    window_size = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % window_size)
    driver = webdriver.Chrome(chrome_options=chrome_options)
    errors_fpath = config["ErrorsFilePath"]
    for gene_name in gene_list:
        try:
            alias_GC_query(driver, gene_name)
        except GeneCardsError as gc_error:
            write_errors(errors_fpath,gene_name,gc_error)


from selenium import webdriver
driver = webdriver.Chrome()
driver.get("https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC4A2")