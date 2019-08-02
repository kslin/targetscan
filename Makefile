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