clean:
	find . -name '*.pyc' -delete


refseq_utr90:
	mkdir -p /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90
	python2 calculate_bins/calculate_bins.py \
		--utr3_file /lab/solexa_bartel/klin/miRNA_models_data/transcript_info/refseq_utr90/UTR_Sequences.txt \
		--tree_file PCT_parameters/Tree.generic.txt \
		--out /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/bins.txt \
		--ref_species 9606 \
		--futures
	python2 find_sites/find_sites.py \
		--seed_file /lab/solexa_bartel/klin/miRNA_models_data/miRNA_info/charlie_mirs/seed_file.txt \
		--bin_file /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/bins.txt \
		--utr3_file /lab/solexa_bartel/klin/miRNA_models_data/transcript_info/refseq_utr90/UTR_Sequences.txt \
		--tree_path PCT_parameters/ \
		--out /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/sites.txt \
		--ribosome_shadow 14 \
		--ref_species 9606 \
		futures
	python2 calculate_features/calculate_features.py \
		--mirna_file /lab/solexa_bartel/klin/miRNA_models_data/miRNA_info/charlie_mirs/mirna_file.txt \
		--site_file /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/sites.txt \
		--ta_sps_file calculate_features/TA_SPS_by_seed_region.txt \
		--orf_file /lab/solexa_bartel/klin/miRNA_models_data/transcript_info/refseq_utr90/orf90.txt \
		--rnaplfold_temp /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/RNAPLFOLD_TEMP \
		--out /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/features.txt \
		--ref_species 9606 \
		--futures
	rm -rf /lab/solexa_bartel/klin/miRNA_models_data/ts7_outputs/refseq_utr90/RNAPLFOLD_TEMP


tests:
	mkdir -p /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated
	python2 calculate_bins/calculate_bins.py \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_file PCT_parameters/Tree.generic.txt \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt \
		--ref_species 9606
	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/bins.txt
	python2 find_sites/find_sites.py \
		--seed_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/seed_file.txt \
		--bin_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/bins.txt \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_path PCT_parameters/ \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt \
		--ribosome_shadow 14 \
		--ref_species 9606
	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/sites.txt
	python2 calculate_features/calculate_features.py \
		--mirna_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/mirna_file.txt \
		--site_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt \
		--ta_sps_file calculate_features/TA_SPS_by_seed_region.txt \
		--orf_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/ORF_Sequences_Human.txt \
		--rnaplfold_temp /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features.txt \
		--ref_species 9606
	rm -rf /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP
	sed 's/[^\t]*\t//28' /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features.txt > /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features_noairs.txt
	diff /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features_noairs.txt /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/original/features.txt

tests_mouse:
	mkdir -p /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse
	python2 calculate_bins/calculate_bins.py \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_file PCT_parameters/Tree.generic.txt \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/bins.txt \
		--ref_species 10090 \
		--futures
	python2 find_sites/find_sites.py \
		--seed_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/seed_file.txt \
		--bin_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/bins.txt \
		--utr3_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/UTR_Sequences.txt \
		--tree_path PCT_parameters/ \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/sites.txt \
		--ribosome_shadow 14 \
		--ref_species 10090 \
		--futures
	python2 calculate_features/calculate_features.py \
		--mirna_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/mirna_file_mouse.txt \
		--site_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/sites.txt \
		--ta_sps_file calculate_features/TA_SPS_by_seed_region.txt \
		--orf_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/ORF_Sequences_Mouse.txt \
		--rnaplfold_temp /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/RNAPLFOLD_TEMP \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/features.txt \
		--ref_species 10090 \
		--futures
	rm -rf /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse/RNAPLFOLD_TEMP

temp:
	python2 calculate_predictions/calculate_predictions.py \
		--feature_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features_with_airs.txt \
		--coeff_file calculate_predictions/Agarwal_2015_parameters.txt \
		--out_dir /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated_mouse



features_with_airs:
	python2 calculate_features/calculate_features.py \
		--mirna_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/mirna_file.txt \
		--site_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/sites.txt \
		--ta_sps_file calculate_features/TA_SPS_by_seed_region.txt \
		--orf_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/ORF_Sequences_Human.txt \
		--rnaplfold_temp /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP \
		--out /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features_with_airs.txt \
		--ref_species 9606 \
		--airs_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/inputs/AIRs_all_cell_lines.txt 
	rm -rf /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/RNAPLFOLD_TEMP

scores:
	mkdir -p /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated2
	python2 calculate_predictions/calculate_predictions.py \
		--feature_file /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated/features_with_airs.txt \
		--coeff_file calculate_predictions/Agarwal_2015_parameters.txt \
		--out_dir /lab/bartel4_ata/kathyl/TargetScan/targetscan/testing/outputs/updated2
