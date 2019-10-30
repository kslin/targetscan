import find_sites_helpers as fsh


def get_all_site_info(gene, group, utr_df, seed_to_species, ref_species, species_to_path, tree_bin_dict, pct_params):
    """
    Given a list of aligned sequences of a single gene, find the site
    information for all the miRNAs with sites in this gene

    Parameters:
    ----------
    gene: string, name of gene

    group: pandas DataFrame, dataframe of miRNA information

    utr_df: pandas DataFrame, dataframe of aligned sequences

    seed_to_species: dictionary, links seed to list of species with that miRNA

    ref_species: string, reference species

    species_to_path: dictionary, tree path for each species

    tree_bin_dict: dictionary, bin to bin-specific tree

    pct_params: DataFrame, PCT parameters

    Output:
    ------
    list of lists, each list has site information for one miRNA
    """

    data = []

    # iterate through miRNAs
    for row in group.iterrows():

        # extract seed, utr, etc information
        seed = row[1]['Seed']
        family = row[1]['miRNA family']
        utr_no_gaps = row[1]['UTR sequence']
        num_sites = row[1]['Num sites']
        bin_specific_tree = tree_bin_dict[row[1]['Bin']]
        species_with_mirna = seed_to_species[seed]
        site_starts, site_ends, site_types = fsh.get_site_info(utr_no_gaps,
                                                               seed)

        # make sure we found the same number of sites as we did before
        assert(num_sites == len(site_starts))

        # find species that have a site at the same place as the ref species
        aligning_species = fsh.find_aligning_species(utr_df, seed,
                                                     species_with_mirna, ref_species,
                                                     num_sites,
                                                     site_types)

        # calculate branch length scores and PCTs
        blss, pct, conserved = zip(*[fsh.calculate_pct(aligning, bin_specific_tree,
                                                       site_type, seed, pct_params)
                                     for (aligning, site_type)
                                     in zip(aligning_species, site_types)])

        # calculate un-normalized branch length scores
        uw_blss = [fsh.get_branch_length_score_generic(aligning, species_to_path)
                   for aligning in aligning_species]

        zipped = zip([gene]*num_sites,
                     [family]*num_sites,
                     [seed]*num_sites,
                     site_starts,
                     site_ends,
                     site_types,
                     pct,
                     conserved,
                     blss,
                     uw_blss,
                     [tuple(x) for x in aligning_species])

        zipped = [list(x) for x in zipped]

        data += zipped

    return data
