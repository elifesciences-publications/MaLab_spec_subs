
def find_ref_seqs(gene_name, tsv_df, errors_fpath):
    """Returns a list of the orthodb ids of the reference sequences from an OrthoDB tsv_df and a set containing
    gene_name and the GeneCards primary alias for gene_name (if it differs from gene_name)
    These reference sequences are defined as records with pub_gene_id, og_name, or description having a
    text match to either gene_name or one of the GeneCards listed aliases for gene_name.
    Alias data is fetched from GeneCards automatically and stored in the aliases directory as text file lists.
    """
    ref_ids = []
    # Go to genecards page for gene_name, extract information for aliases from the webpage/ html
    aliases_fpath = "aliases_data/" + gene_name + "_aliases.txt"

    # If file doesn't exist or was improperly downloaded to yield only one line, repeat
    # fetching alias names
    if (not os.path.exists(aliases_fpath)) or len(open(aliases_fpath, 'r').readlines()) == 1:
        try:
            aliases, gc_name = alias_GC_query(gene_name)
            matches = set((gene_name.upper(), gc_name.upper()))
        except GeneCardsError as gc_error:
            aliases = [gene_name]
            matches = set((gene_name.upper(),))
            write_errors(errors_fpath, gene_name, gc_error)
    else:
        # Read aliases information previously downloaded from GeneCards
        aliases_f = open(aliases_fpath, 'r')
        aliases = aliases_f.readlines()
        aliases_f = open(aliases_fpath, 'r')
        gc_name = aliases_f.readline().strip()
        matches = set((gene_name.upper(), gc_name.upper()))
    # Remove spaces, commas, new line chars, and capitalization from alias strings
    formatted_aliases = [format_odb_field(alias) for alias in aliases]
    # Search fields in search_fields for matches to the alias strings provided by GeneCards
    # Iterate tsv_df rows, save all reference ids which have matches
    search_fields = ["pub_gene_id", "og_name", "description"]
    aliases_pat = "|".join(formatted_aliases)
    for idx, row in tsv_df.iterrows():
        #         for field in search_fields:
        # Current behavior: exact matches in formatted pub_gene_id, og_name, or description only.
        # TODO: Add in partial string matching. Difficulties with distinguishing gene_names
        for alias in formatted_aliases:
            for field in search_fields:
                formatted_field = format_odb_field(str(row[field]))
                try:
                    if re.search(alias, formatted_field):
                        if idx not in ref_ids:
                            ref_ids.append(idx)  # ["int_prot_id"])
                            break
                except Exception as e:
                    # Special regexp characters present in alias
                    if alias in formatted_field:
                        if idx not in ref_ids:
                            ref_ids.append(idx)  # ["int_prot_id"])
                            break
    # matches is a tuple of strings used to filter reference sequences in pg_id_df; either
    # one entry or gene_name and then the entry used by the GeneCards page
    return ref_ids, matches
