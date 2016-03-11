clean:
	find . -name '*.pyc' -delete


bins:
	# cd calculate_bins && python calculate_bins.py ../../infiles/UTR_Sequences_sample.txt ../../outfiles/bins_tiny_outfile.txt
	cd calculate_bins && python calculate_bins.py ../../infiles/tinier_utrs.txt ../../outfiles/bins_tinier_outfile.txt
	# cd calculate_bins && python calculate_bins.py ../../infiles/tiny_utrs.txt ../../outfiles/bins_small_outfile.txt
	# cd calculate_bins && python calculate_bins.py ../../infiles/UTR_Sequences_Ensembl_forTS_kathy_no_header.txt ../../outfiles/bins_big_outfile.txt

sites:
	# cd find_sites && python find_sites.py ../../infiles/miR_Family_info_sample.txt ../../outfiles/bins_tiny_outfile.txt ../../infiles/UTR_Sequences_sample.txt ../../outfiles/sites_tiny_outfile.txt
	cd find_sites && python find_sites.py ../../infiles/kathy_seed_file.txt ../../outfiles/bins_tinier_outfile.txt ../../infiles/tinier_utrs.txt ../../outfiles/sites_tinier_outfile.txt
	# cd find_sites && python find_sites.py ../../infiles/kathy_seed_file.txt ../../outfiles/bins_small_outfile.txt ../../infiles/tiny_utrs.txt ../../outfiles/sites_small_outfile.txt
	# cd find_sites && python find_sites.py ../../infiles/kathy_seed_file.txt ../../outfiles/bins_big_outfile.txt ../../infiles/UTR_Sequences_Ensembl_forTS_kathy_no_header.txt ../../outfiles/sites_big_outfile.txt

features:
	# cd calculate_features/ && python calculate_features.py  ../../infiles/miR_for_context_scores.sample.txt ../../outfiles/sites_tiny_outfile.txt ../../infiles/ORF_Sequences_sample.txt ../../outfiles/tiny_output.txt
	cd calculate_features/ && python calculate_features.py ../../infiles/kathy_mirna_file.txt ../../outfiles/sites_tinier_outfile.txt ../../infiles/ORF_Sequences_human.txt ../../outfiles/tinier_output.txt
	# cd calculate_features/ && python calculate_features.py ../../infiles/kathy_mirna_file.txt ../../outfiles/sites_small_outfile.txt ../../infiles/ORF_Sequences_human.txt ../../outfiles/small_output.txt
	# cd calculate_features/ && python calculate_features.py ../../infiles/kathy_mirna_file.txt ../../outfiles/sites_big_outfile.txt ../../infiles/ORF_Sequences_human.txt ../../outfiles/big_output.txt