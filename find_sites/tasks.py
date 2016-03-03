import find_sites_helpers

# def get_all_site_info(gene, family, seed, utr_df, species, num_sites, utr_no_gaps, window_dict):
#     site_starts, site_ends, site_types = find_sites_helpers.get_site_info(utr_no_gaps, seed, window_dict)
#     assert(num_sites == len(site_starts))
#     aligning_species = find_sites_helpers.find_aligning_species(utr_df,seed,species,num_sites)

#     zipped = zip([gene]*num_sites, [family]*num_sites, [seed]*num_sites, site_starts, site_ends, site_types, aligning_species)
#     zipped = [list(x) for x in zipped]

#     return zipped

def get_all_site_info(gene,group,utr_df,SEED_TO_SPECIES,SEED_WINDOW_DICT):

    data = []
    for row in group.iterrows():
        seed = row[1]['Seed']
        family = row[1]['miRNA family']
        utr_no_gaps = row[1]['UTR sequence']
        num_sites = row[1]['Num sites']
        # bin = row[1]['Bin']
        window_dict = SEED_WINDOW_DICT[seed]
        species = SEED_TO_SPECIES[seed]
        site_starts, site_ends, site_types = find_sites_helpers.get_site_info(utr_no_gaps, seed, window_dict)
        assert(num_sites == len(site_starts))
        aligning_species = find_sites_helpers.find_aligning_species(utr_df,seed,species,num_sites)

        # pct,conserved = zip(*[find_sites_helpers.calculate_pct(aligning,bin,site_type,seed) for (aligning,site_type) in zip(aligning_species,site_types)])

        # zipped = zip([gene]*num_sites, [family]*num_sites, [seed]*num_sites, site_starts, site_ends, site_types, pct, conserved)
        # zipped = [list(x) for x in zipped]

        blss = [find_sites_helpers.get_branch_length_score(aligning) for aligning in aligning_species]

        zipped = zip([gene]*num_sites, [family]*num_sites, [seed]*num_sites, site_starts, site_ends, site_types, blss)
        zipped = [list(x) for x in zipped]

        data += zipped

    return data