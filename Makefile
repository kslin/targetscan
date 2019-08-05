clean:
	find . -name '*.pyc' -delete


bins:
	mkdir -p /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated
	python2 calculate_bins/calculate_bins.py \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_file PCT_parameters/Tree.generic.txt \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt \
		--ref_species 9606 \
		--futures
	sort "-k1,1" /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt -o /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt
	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/bins.txt


sites:
	python2 find_sites/find_sites.py \
		--seed_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/seed_file.txt \
		--bin_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_path PCT_parameters/ \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt \
		--ribosome_shadow 14 \
		--ref_species 9606
	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/sites.txt


features:
	python2 calculate_features/calculate_features.py \
		--mirna_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/mirna_file.txt \
		--site_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt \
		--ta_sps_file calculate_features/TA_SPS_by_seed_region.txt \
		--orf_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/ORF_Sequences_Human.txt \
		--rnaplfold_temp /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features.txt \
		--ref_species 9606 \
		--futures
	rm -rf /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP
# 	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/features.txt