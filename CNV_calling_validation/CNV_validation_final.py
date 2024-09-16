import numpy as np
from collections import defaultdict
import sys

strain = sys.argv[1]

FILTERED_BAM_DIR = '/storage/NFS/ANALYSIS/DNAseq/PfResistome_2023/bam_sam_filtered_daisy_%s' % strain
OUTPUT_DIR = '/storage/NFS/ANALYSIS/DNAseq/PfResistome_2023/daisy_test' # '/projects/winzeler/ANALYSIS/DNAseq/test/daisy_test'

WINDOW = 10000

# ================================================================
# Resistome dataset metadata
# ================================================================

bad_good_clone_dict = {'11C-R-385-GCCAAT': '11C-R-385',
'11C-R-474-ACTTGA': '11C-R-474',
'6A-R-385-ACAGTG': '6A-R-385',
'6A-R-474-CAGATC': '6A-R-474',
'Wirth-TCMDC125334-3D7-S2-F1-A10': 'TCMDC125334-3D7-T-F1-C1',
'Wirth-TCMDC125334-3D7-S2-F1-B6': 'TCMDC125334-3D7-T-F1-C2',
'Wirth-TCMDC125334-3D7-S2-F1-C7': 'TCMDC125334-3D7-T-F1-C3',
'Wirth-TCMDC125334-3D7-S2-F2-A3': 'TCMDC125334-3D7-T-F2-C1',
'Wirth-TCMDC125334-3D7-S2-F2-A6': 'TCMDC125334-3D7-T-F2-C2',
'Wirth-TCMDC125334-3D7-S2-F2-B9': 'TCMDC125334-3D7-T-F2-C3',
'Wirth-TCMDC125334-3D7-S2-F2-C10': 'TCMDC125334-3D7-T-F2-C4',
'Wirth-TCMDC125334-3D7-S2-F2-C8': 'TCMDC125334-3D7-T-F2-C5'}

def bam_name_to_clone_name(sample_bam_name):
	sample_name = sample_bam_name.split('_p')[0].split('_s')[0]
	if sample_name.startswith('p_fal_'):
		sample_name = sample_name.split('p_fal_')[1]
	
	if 'Actinomycin' in sample_name:
		sample_name = sample_name.replace('Actinomycin', 'Antimycin')
	
	if sample_name in bad_good_clone_dict:
		sample_name = bad_good_clone_dict[sample_name]
	
	if 'NCP26' in sample_name and 'clone' in sample_name:
		flask = sample_name.split('-F')[1][0]
		well = sample_name.split('-clone')[1]
		sample_name = 'Wirth-NCP26-Dd2-' + flask + well
	
	return sample_name

# ================================================================
# Load core vs. non-core regions
# ================================================================

bangladesh_dir = '/storage/NFS/ROTATION_PROJECT/daisy/bangladesh'
f = open('%s/bangladesh_samples.txt' % bangladesh_dir, 'r')
bangladesh_samples = [line.strip() for line in f]
vcf_file = open('%s/data/%s.ann.txt' % (bangladesh_dir, bangladesh_samples[0]), 'rt')

chrom_length_dict = {}

for line in vcf_file:
    if line[:8] == '##contig':
        chrom = line.split('ID=')[1].split(',')[0]
        length = int(line.split('length=')[1].split('>')[0])
        chrom_length_dict[chrom] = length
    if line[:2] != '##':
        break

chromosomes = list(chrom_length_dict.keys())

genome_ref_dir = '/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA/pfal/'
f = open("%s/3D7_chromosome_noncore_regions.txt" % (genome_ref_dir), 'r')

chrom_noncore_regions_dict = defaultdict(list)

for line in f:
	items = line.strip().split(' ')
	chrom, start, end = items
	if chrom not in ['M76611', 'PFC10_API_IRAB']:
		if start == 'start':
			start = 1
		elif end == 'end':
			end = chrom_length_dict[chrom]
		
		end = int(end); start = int(start)
		chrom_noncore_regions_dict[chrom].append((start, end))

NONCORE_BUFFER = 10000

def pos_near_noncore_region(chrom, pos):
	for start, end in chrom_noncore_regions_dict[chrom]:
		if pos > (start - NONCORE_BUFFER) and pos < (end + NONCORE_BUFFER):
			return (start, end)
	return False

def pos_near_subtelomere(chrom, pos, buffer = NONCORE_BUFFER):
	for start, end in chrom_noncore_regions_dict[chrom][:2]:
		if pos > (start - buffer) and pos < (end + buffer):
			return True
	return False

# ================================================================
# Get copy ratio data
# ================================================================

def get_copyratio_dict(sample_full_name, strain):
	chrom_interval_copyratio_dict = defaultdict(dict)
	
	copyratio_dir = '/storage/NFS/ANALYSIS/DNAseq/PfResistome_2023/gatk_cnv_%s' % strain
	f = open('%s/%s.denoisedCR.tsv' % (copyratio_dir, sample_full_name), 'r')
	line = f.readline()
	while line[0] == '@':
		line = f.readline()
	
	header = line.strip().split('\t')
	for line in f:
		chrom, start, end, log2_copyratio = line.strip().split('\t')
		chrom_interval_copyratio_dict[chrom][(int(start), int(end))] = float(log2_copyratio)
	
	return chrom_interval_copyratio_dict

def get_copyratio_list(chrom_interval_copyratio_dict, chrom, start_pos, end_pos):
	copyratio_list = []
	start_flag = False
	for interval_start, interval_end in sorted(chrom_interval_copyratio_dict[chrom]):
		copyratio = chrom_interval_copyratio_dict[chrom][(interval_start, interval_end)]
		if start_pos >= interval_start and start_pos <= interval_end:
			copyratio_list.append(copyratio)
			start_flag = True
		if end_pos >= interval_start and end_pos <= interval_end:
			copyratio_list.append(copyratio)
			break
		if start_flag:
			copyratio_list.append(copyratio)
	return copyratio_list

# ================================================================
# Get coverage data
# ================================================================

def get_sample_coverage_dicts(strain):
	sample_coverage_dict = defaultdict(int)
	sample_avg_coverage_dict = defaultdict(float)
	sample_perc_bases_gt_20_dict = defaultdict(float)
	
	output_dir = '/storage/NFS/ANALYSIS/DNAseq/PfResistome_2023/daisy_test'
	f = open('%s/coverage_summary_%s.txt' % (output_dir, strain), 'r')
	header = f.readline()
	for line in f:
		sample_bam_name, sample_name, total_cov, avg_cov, num_bases_gt_20, perc_bases_gt_20 = line.strip().split('\t')
		sample_coverage_dict[sample_name] = int(total_cov)
		sample_perc_bases_gt_20_dict[sample_name] = float(perc_bases_gt_20)
		sample_avg_coverage_dict[sample_name] = float(avg_cov)
	
	return sample_coverage_dict, sample_avg_coverage_dict, sample_perc_bases_gt_20_dict

# ================================================================
# Get predicted CNVs
# ================================================================

def get_predicted_CNVs(sample_name):
	predicted_CNVs = []
	
	f = open('%s/NewSummmaryTableall_with_clado.tsv' % OUTPUT_DIR, 'r')
	header = f.readline()
	for line in f:
		items = line.rstrip().split('\t')
		if items[0] == "" or items[1] == "None found":
			continue
		clone_name, chrom, CNV_size, avg_val, CNV_start, CNV_end, CNV_num, pval, IQR, parent_avg, strain = items
		pval = float(pval)
		if clone_name == sample_name:
			predicted_CNVs.append((chrom, int(CNV_start), int(CNV_end)))
	
	return predicted_CNVs

def get_clone_compound_strain_dict():
	clone_compound_dict = {}; clone_strain_dict = {}
	f = open('%s/NewSummmaryTableall_with_clado.tsv' % OUTPUT_DIR, 'r')
	header = f.readline()
	for line in f:
		items = line.rstrip().split('\t')
		if items[0] == "" or items[1] == "None found":
			continue
		clone_name, chrom, CNV_size, avg_val, CNV_start, CNV_end, CNV_num, pval, IQR, parent_avg, strain = items
		clone_compound_dict[clone_name] = 'N/A'
		clone_strain_dict[clone_name] = strain
	
	return clone_compound_dict, clone_strain_dict

def get_CNV_pval_dict():
	CNV_pval_dict = {}	
	f = open('%s/NewSummmaryTableall_with_clado.tsv' % OUTPUT_DIR, 'r')
	header = f.readline()
	for line in f:
		items = line.rstrip().split('\t')
		if items[0] == "" or items[1] == "None found":
			continue
		clone_name, chrom, CNV_size, avg_val, CNV_start, CNV_end, CNV_num, pval, IQR, parent_avg, strain = items
		CNV_pval_dict[(clone_name, chrom, int(CNV_start), int(CNV_end))] = float(pval)
	
	return CNV_pval_dict

# ================================================================
# Get reads that may serve as evidence for tandem duplication
# ================================================================

def get_supporting_reads_dict(sample_bam_name):
	chrom_reads_dict = {rtype: defaultdict(list) for rtype in ['RR', 'RF', 'FF']}
	
	try:
		f = open('%s/%s_INVDUP.sam' % (FILTERED_BAM_DIR, sample_bam_name), 'r')
	except:
		return False
	
	for line in f:
		items = line.strip().split('\t')
		qname, flag, chrom, pos, mapq, cigar, chrom_next, pos_next, tlen, seq, qual = items[:11]
		flag = int(flag); pos = int(pos); mapq = int(mapq); pos_next = int(pos_next); insert_size = np.abs(int(tlen))
		
		if mapq >= 20:
			if flag in [97, 81, 145, 161] and insert_size > 1000:
				chrom_reads_dict['RF'][chrom].append((flag, pos, len(seq), chrom_next, pos_next))
			elif flag in [177, 113]:
				chrom_reads_dict['RR'][chrom].append((flag, pos, len(seq), chrom_next, pos_next))
			elif flag in [65, 129]:
				chrom_reads_dict['FF'][chrom].append((flag, pos, len(seq), chrom_next, pos_next))
	
	return chrom_reads_dict

# ==================================================================
# Sliding window to find most likely location of start and end pos
# ==================================================================

def infer_tandem_boundaries(start_pos_count_dict, end_pos_count_dict):
	window = 100
	window_start_count_dict = defaultdict(int) # Start of window -> count of read start positions in that window
	
	for ipos in start_pos_count_dict:
		for pos in np.arange(ipos, ipos + window):
			if pos in start_pos_count_dict:
				window_start_count_dict[ipos] += start_pos_count_dict[pos]
	
	best_window_start = max(start_pos_count_dict.keys())
	max_count = window_start_count_dict[best_window_start]
	
	for window_start in sorted(window_start_count_dict.keys(), reverse=True):
		if window_start_count_dict[window_start] > max_count:
			best_window_start = window_start
			max_count = window_start_count_dict[window_start]
	
	start_support = max_count
	mode_pos = best_window_start
	mode = start_pos_count_dict[mode_pos] # Mode within the window
	
	for pos in np.arange(best_window_start, best_window_start + window):
		if start_pos_count_dict[pos] > mode:
			mode_pos = pos
			mode = start_pos_count_dict[pos]
	
	# Starting from mode of window, move left until read count is <=1 within window of 20 on the left side	
	predicted_start_pos = mode_pos # Predicted start of CNV
	cur_read_count = sum([start_pos_count_dict[pos] for pos in range(predicted_start_pos-20, predicted_start_pos)])	
	while cur_read_count > 1:
		predicted_start_pos -= 1
		cur_read_count = sum([start_pos_count_dict[pos] for pos in range(predicted_start_pos-20, predicted_start_pos)])
	
	window = 100
	window_start_count_dict = defaultdict(int) # Start of window -> count of read end positions in that window
	
	for ipos in end_pos_count_dict:
		for pos in np.arange(ipos, ipos + window):
			if pos in end_pos_count_dict:
				window_start_count_dict[ipos] += end_pos_count_dict[pos]
	
	best_window_start = min(end_pos_count_dict.keys())
	max_count = window_start_count_dict[best_window_start]
	
	for window_start in window_start_count_dict:
		if window_start_count_dict[window_start] > max_count:
			best_window_start = window_start
			max_count = window_start_count_dict[window_start]
	
	end_support = max_count
	mode_pos = best_window_start
	mode = end_pos_count_dict[mode_pos] # Mode within the window
	
	for pos in np.arange(best_window_start, best_window_start + window):
		if end_pos_count_dict[pos] > mode:
			mode_pos = pos
			mode = end_pos_count_dict[pos]
	
	# Starting from mode of window, move right until read count is <=1 within window of 20 on the right side
	predicted_end_pos = mode_pos
	cur_read_count = sum([end_pos_count_dict[pos] for pos in range(predicted_end_pos+1, predicted_end_pos+21)])	
	while cur_read_count > 1:
		predicted_end_pos += 1
		cur_read_count = sum([end_pos_count_dict[pos] for pos in range(predicted_end_pos+1, predicted_end_pos+21)])
	
	return (predicted_start_pos, predicted_end_pos, start_support, end_support)

# ==================================================================
# Verify consistent mapping of mates
# ==================================================================

def check_consistent_mates(pos_pair_dict, filter_func=lambda x: True, threshold=5, window=50, single_chrom=False):
		
	chrom_pos_count_dict = {chrom: defaultdict(int) for chrom in chromosomes}
	subtelomere_count = 0
	for pos in pos_pair_dict:
		if filter_func(pos):
			for key in pos_pair_dict[pos]:
				if single_chrom:
					mate_pos = key; mate_chrom = single_chrom
				else:
					mate_chrom, mate_pos = key
				if pos_near_subtelomere(mate_chrom, mate_pos, buffer = 100):
					subtelomere_count += 1
				else:
					chrom_pos_count_dict[mate_chrom][mate_pos] += 1
	
	if subtelomere_count >= threshold:
		return ("Subtelomere", -1, subtelomere_count)
	
	for chrom in chrom_pos_count_dict:
		window_start_count_dict = defaultdict(int) # Start of window -> count of read positions in that window
		for ipos in chrom_pos_count_dict[chrom].keys():
			for pos in np.arange(ipos, ipos + window):
				if pos in chrom_pos_count_dict[chrom]:
					window_start_count_dict[ipos] += chrom_pos_count_dict[chrom][pos]
		if len(window_start_count_dict) == 0:
			continue
		max_count_pos, max_count = sorted(window_start_count_dict.items(), key=lambda x: x[1])[-1]
		if max(window_start_count_dict.values()) >= threshold:
			return (chrom, max_count_pos, max_count)
	
	return False

# ==================================================================
# Sliding window to find potential enrichment peak
# ==================================================================

def quantify_enrichment(chrom_pos_count_dict, filter_func, window=50):
	if len(chrom_pos_count_dict) == 0:
		return (0, 0, 0)
	
	window_start_count_dict = defaultdict(int) # Start of window -> count of read positions in that window

	for ipos in chrom_pos_count_dict.keys():
		if filter_func(ipos):
			for pos in np.arange(ipos, ipos + window):
				if pos in chrom_pos_count_dict:
					window_start_count_dict[ipos] += chrom_pos_count_dict[pos]

	best_window_start = min(chrom_pos_count_dict.keys())
	max_count = window_start_count_dict[best_window_start]

	for window_start in window_start_count_dict:
		if window_start_count_dict[window_start] > max_count:
			best_window_start = window_start
			max_count = window_start_count_dict[window_start]

	mode_pos = best_window_start
	mode = chrom_pos_count_dict[mode_pos] # Mode within the window
		
	for pos in np.arange(best_window_start, best_window_start + window):
		if chrom_pos_count_dict[pos] > mode:
			mode_pos = pos
			mode = chrom_pos_count_dict[pos]
	
	return (max_count, mode, mode_pos)

# ================================================================
# Filter pos_count_dict to positions within window with the highest total count
# ================================================================

def filter_pos_count_dict_to_window(pos_count_dict, window=300):
	
	if len(pos_count_dict) == 0:
		return pos_count_dict
	
	sorted_pos_list = sorted(pos_count_dict.keys())
	window_count_dict = defaultdict(int)
	
	for i in range(len(sorted_pos_list)):
		start_pos = sorted_pos_list[i]; j = i
		while j < len(sorted_pos_list) and sorted_pos_list[j] < start_pos + window:
			window_count_dict[start_pos] += pos_count_dict[sorted_pos_list[j]]; j += 1
	
	max_count_start_pos, max_count = sorted(window_count_dict.items(), key=lambda x: x[1])[-1]
	
	new_pos_count_dict = defaultdict(int)
	for pos in sorted_pos_list:
		if pos >= max_count_start_pos and pos < (max_count_start_pos + window):
			new_pos_count_dict[pos] = pos_count_dict[pos]
	
	return new_pos_count_dict

# ================================================================
# Validate a particular CNV
# ================================================================

def pos_within_window(pos, ref_pos, window):
	return pos >= (ref_pos - window) and pos <= (ref_pos + window)

def validate_CNV(cnv_info, chrom_reads_dict):
	
	RF_start_pos_count_dict = defaultdict(int)
	RF_end_pos_count_dict = defaultdict(int)
	
	RR_pos_pair_dict = defaultdict(list) # Pos within window of either CNV start or end, mate pos
	RR_pos_count_dict = defaultdict(int)
	FF_pos_pair_dict = defaultdict(list) # Pos within window of either CNV start or end, mate pos
	FF_pos_count_dict = defaultdict(int)
	
	RR_otherchrom_pos_pair_dict = defaultdict(list) # Pos within window of either CNV start or end, (mate chrom, mate pos)
	RR_otherchrom_pos_count_dict = defaultdict(int)
	FF_otherchrom_pos_pair_dict = defaultdict(list) # Pos within window of either CNV start or end, (mate chrom, mate pos)
	FF_otherchrom_pos_count_dict = defaultdict(int)
	
	tandem_supporting_reads = 0
	total_supporting_reads = 0
	
	tandem_altered_CNV_boundary_flag = False
	
	CNV_chrom, CNV_start, CNV_end = cnv_info
	
	pos_in_start_candidate_region = lambda pos: pos_within_window(pos, CNV_start, WINDOW)
	pos_in_end_candidate_region = lambda pos: pos_within_window(pos, CNV_end, WINDOW)
	
	# Check if either start or end is close to a hypervariable region
	start_noncore_region = pos_near_noncore_region(CNV_chrom, CNV_start)
	end_noncore_region = pos_near_noncore_region(CNV_chrom, CNV_end)
	
	if start_noncore_region is not False: # Start near noncore
		noncore_start, noncore_end = start_noncore_region
		pos_in_start_candidate_region = lambda pos: (pos >= noncore_start - 1000 and pos <= noncore_end + 1000) or pos_within_window(pos, CNV_start, WINDOW)
		
	if end_noncore_region is not False: # End near noncore
		noncore_start, noncore_end = end_noncore_region
		pos_in_end_candidate_region = lambda pos: (pos >= noncore_start - 1000 and pos <= noncore_end + 1000) or pos_within_window(pos, CNV_end, WINDOW)
	
	# Validate each CNV: start by searching for potentially supporting reads
	for flag, pos, seq_len, mate_chrom, mate_pos in chrom_reads_dict['RF'][CNV_chrom]:
		
		if flag in [81, 145] and pos_in_start_candidate_region(pos):
			if mate_chrom == '=' and pos_in_end_candidate_region(mate_pos + seq_len):
				RF_start_pos_count_dict[pos] += 1
				RF_end_pos_count_dict[mate_pos + seq_len] += 1
				tandem_supporting_reads += 1
		
		elif flag in [97, 161] and pos_in_end_candidate_region(pos + seq_len):
			if mate_chrom == '=' and pos_in_start_candidate_region(mate_pos):
				RF_start_pos_count_dict[mate_pos] += 1
				RF_end_pos_count_dict[pos + seq_len] += 1
				tandem_supporting_reads += 1
	
	if tandem_supporting_reads <= 2: # Only if there are not enough initial candidates, consider incorrect CNV segmentation
		tandem_altered_CNV_boundary_flag = True
		RF_start_pos_count_dict = defaultdict(int)
		RF_end_pos_count_dict = defaultdict(int)
		for flag, pos, seq_len, mate_chrom, mate_pos in chrom_reads_dict['RF'][CNV_chrom]:
			if flag in [81, 145] and pos_within_window(pos, CNV_start, WINDOW) and mate_chrom == '=' and (mate_pos - pos) > 1000:
					RF_start_pos_count_dict[pos] += 1
					RF_end_pos_count_dict[mate_pos + seq_len] += 1
		RF_start_pos_count_dict = filter_pos_count_dict_to_window(RF_start_pos_count_dict, window=300)
		RF_end_pos_count_dict = filter_pos_count_dict_to_window(RF_end_pos_count_dict, window=300)
		tandem_supporting_reads = sum(RF_start_pos_count_dict.values()) # Override
	
	if tandem_supporting_reads <= 2: # Only if there are not enough initial candidates, consider incorrect CNV segmentation
		tandem_altered_CNV_boundary_flag = True
		RF_start_pos_count_dict = defaultdict(int)
		RF_end_pos_count_dict = defaultdict(int)
		for flag, pos, seq_len, mate_chrom, mate_pos in chrom_reads_dict['RF'][CNV_chrom]:
			if flag in [97, 161] and pos_within_window(pos + seq_len, CNV_end, WINDOW) and mate_chrom == '=' and (pos - mate_pos) > 1000:
				RF_start_pos_count_dict[mate_pos] += 1
				RF_end_pos_count_dict[pos + seq_len] += 1
		RF_start_pos_count_dict = filter_pos_count_dict_to_window(RF_start_pos_count_dict, window=300)
		RF_end_pos_count_dict = filter_pos_count_dict_to_window(RF_end_pos_count_dict, window=300)
		tandem_supporting_reads = sum(RF_start_pos_count_dict.values()) # Override
	
	total_supporting_reads += sum(RF_start_pos_count_dict.values())
	
	for flag, pos, seq_len, mate_chrom, mate_pos in chrom_reads_dict['RR'][CNV_chrom]:
		if (pos_within_window(pos, CNV_start, WINDOW) or pos_within_window(pos, CNV_end, WINDOW)):
			if mate_chrom == '=':
				RR_pos_pair_dict[pos].append(mate_pos)
				RR_pos_count_dict[pos] += 1
			else:
				RR_otherchrom_pos_pair_dict[pos].append((mate_chrom, mate_pos))
				RR_otherchrom_pos_count_dict[pos] += 1
			total_supporting_reads += 1
		
		if mate_chrom == '=' and ((pos_within_window(mate_pos, CNV_start, WINDOW) or pos_within_window(mate_pos, CNV_end, WINDOW))):
			RR_pos_pair_dict[mate_pos].append(pos)
			RR_pos_count_dict[mate_pos] += 1
			total_supporting_reads += 1
	
	for flag, pos, seq_len, mate_chrom, mate_pos in chrom_reads_dict['FF'][CNV_chrom]:
		if (pos_within_window(pos, CNV_start, WINDOW) or pos_within_window(pos, CNV_end, WINDOW)):
			if mate_chrom == '=':
				FF_pos_pair_dict[pos].append(mate_pos)
				FF_pos_count_dict[pos] += 1
			else:
				FF_otherchrom_pos_pair_dict[pos].append((mate_chrom, mate_pos))
				FF_otherchrom_pos_count_dict[pos] += 1
			total_supporting_reads += 1
		
		if mate_chrom == '=' and ((pos_within_window(mate_pos, CNV_start, WINDOW) or pos_within_window(mate_pos, CNV_end, WINDOW))):
			FF_pos_pair_dict[mate_pos].append(pos)
			FF_pos_count_dict[mate_pos] += 1
			total_supporting_reads += 1
	
	# Validation failed	
	if total_supporting_reads < 6 and not tandem_supporting_reads >= 2:
		return ('LACK OF READS', None, None, None)
	
	# ===========================================================
	# Infer type of CNV
	# ===========================================================
	
	if tandem_supporting_reads >= 4: # Potential tandem duplication
		predicted_start_pos, predicted_end_pos, start_support, end_support = infer_tandem_boundaries(RF_start_pos_count_dict, RF_end_pos_count_dict)
		if tandem_altered_CNV_boundary_flag: # Limit to stricter support requirement
			if start_support >= 5 and end_support >= 5:
				return ('TANDEM DUPLICATION', predicted_start_pos, predicted_end_pos, (start_support, end_support))
			else:
				return ('DUBIOUS TANDEM DUPLICATION', predicted_start_pos, predicted_end_pos, (start_support, end_support))
		else:
			return ('TANDEM DUPLICATION', predicted_start_pos, predicted_end_pos, (start_support, end_support))
	
	elif pos_near_subtelomere(CNV_chrom, CNV_start, buffer=75000): # Potentially goes to end of chromosome
		consistent_mates_result = check_consistent_mates(FF_otherchrom_pos_pair_dict, lambda pos: pos_within_window(pos, CNV_end, WINDOW), threshold=6, window=50)
		if consistent_mates_result is not False:
			mate_chrom, mate_max_count_pos, mate_max_count = consistent_mates_result
			max_count, mode, mode_pos = quantify_enrichment(FF_otherchrom_pos_count_dict, lambda pos: pos_within_window(pos, CNV_end, WINDOW), window=50)
			if max_count >= 5 or mode >= 3:
				return ('INTERCHROMOSOMAL, INSERTED AT %s:%i' % (mate_chrom, mate_max_count_pos), 1, mode_pos, max_count)
	
	elif pos_near_subtelomere(CNV_chrom, CNV_end, buffer=75000): # Potentially goes to end of chromosome
		consistent_mates_result = check_consistent_mates(RR_otherchrom_pos_pair_dict, lambda pos: pos_within_window(pos, CNV_start, WINDOW), threshold=6, window=50)
		if consistent_mates_result is not False:
			mate_chrom, mate_max_count_pos, mate_max_count = consistent_mates_result
			max_count, mode, mode_pos = quantify_enrichment(RR_otherchrom_pos_count_dict, lambda pos: pos_within_window(pos, CNV_start, WINDOW), window=50)
			if max_count >= 5 or mode >= 3:
				return ('INTERCHROMOSOMAL, INSERTED AT %s:%i' % (mate_chrom, mate_max_count_pos), mode_pos, chrom_length_dict[CNV_chrom], max_count)
	
	else: # Potentially inverted
		
		consistent_mates_result = check_consistent_mates(RR_otherchrom_pos_pair_dict, lambda pos: pos_in_start_candidate_region(pos), threshold=5, window=50)
		if consistent_mates_result is not False:
			mate_chrom, mate_max_count_pos, mate_max_count = consistent_mates_result
			max_count, mode, mode_pos = quantify_enrichment(RR_otherchrom_pos_count_dict, lambda pos: pos_in_start_candidate_region(pos), window=50)
			if max_count >= 5:
				return ('INTERCHROMOSOMAL INVERTED DUPLICATION, INSERTED AT %s:%i' % (mate_chrom, mate_max_count_pos), mode_pos, None, max_count)
		
		consistent_mates_result = check_consistent_mates(FF_otherchrom_pos_pair_dict, threshold=5, window=50)
		if consistent_mates_result is not False:
			mate_chrom, mate_max_count_pos, mate_max_count = consistent_mates_result
			max_count, mode, mode_pos = quantify_enrichment(RR_otherchrom_pos_count_dict, lambda pos: pos_in_end_candidate_region(pos), window=50)
			if max_count >= 5:
				return ('INTERCHROMOSOMAL INVERTED DUPLICATION, INSERTED AT %s:%i' % (mate_chrom, mate_max_count_pos), None, mode_pos, max_count)
		
		consistent_mates_result = check_consistent_mates(RR_pos_pair_dict, lambda pos: pos_in_start_candidate_region(pos), threshold=5, window=50, single_chrom=CNV_chrom)
		if consistent_mates_result is not False:
			max_count, mode, mode_pos = quantify_enrichment(RR_pos_count_dict, lambda pos: pos_in_start_candidate_region(pos), window=50)
			if max_count >= 8:
				return ('POSSIBLE INVERTED DUPLICATION', mode_pos, None, max_count)
		
		consistent_mates_result = check_consistent_mates(RR_pos_pair_dict, lambda pos: pos_in_end_candidate_region(pos), threshold=5, window=50, single_chrom=CNV_chrom)
		if consistent_mates_result is not False:
			max_count, mode, mode_pos = quantify_enrichment(RR_pos_count_dict, lambda pos: pos_in_end_candidate_region(pos), window=50)
			if max_count >= 8:
				return ('POSSIBLE INVERTED DUPLICATION', None, mode_pos, max_count)
		
		consistent_mates_result = check_consistent_mates(FF_pos_pair_dict, lambda pos: pos_in_start_candidate_region(pos), threshold=5, window=50, single_chrom=CNV_chrom)
		if consistent_mates_result is not False:
			max_count, mode, mode_pos = quantify_enrichment(FF_pos_count_dict, lambda pos: pos_in_start_candidate_region(pos), window=50)
			if max_count >= 8:
				return ('POSSIBLE INVERTED DUPLICATION', mode_pos, None, max_count)
		
		consistent_mates_result = check_consistent_mates(FF_pos_pair_dict, lambda pos: pos_in_end_candidate_region(pos), threshold=5, window=50, single_chrom=CNV_chrom)
		if consistent_mates_result is not False:
			max_count, mode, mode_pos = quantify_enrichment(FF_pos_count_dict, lambda pos: pos_in_end_candidate_region(pos), window=50)
			if max_count >= 8:
				return ('POSSIBLE INVERTED DUPLICATION', None, mode_pos, max_count)
	
	return ("INSUFFICIENT SUPPORT", None, None, -1)

# ================================================================
# Main
# ================================================================

sample_coverage_dict, sample_avg_coverage_dict, sample_perc_bases_gt_20_dict = get_sample_coverage_dicts(strain)
clone_compound_dict, clone_strain_dict = get_clone_compound_strain_dict()
CNV_pval_dict = get_CNV_pval_dict()

output_f = open('%s/CNV_validation_results_clado_%s.tsv' % (OUTPUT_DIR, strain), 'w')
output_f.write('\t'.join(['compound', 'clone_name', 'chrom', 'orig_start', 'orig_end', 'p-value', 'avg_copyratio', 'classification', 'support', 'new_start', 'new_end', 'avg_coverage_genome', 'perc_bases_gt_20_genome', 'prop_RF_genome', 'prop_RRFF_genome']) + '\n')

for clone in ['p_fal_cladorA_p', 'p_fal_cladorB_p', 'p_fal_cladorC_p']:
	clone_compound_dict[clone] = 'cladosporin'
	clone_strain_dict[clone] = 'Dd2'

f = open('/storage/NFS/ANALYSIS/DNAseq/PfResistome_2023/bamnames_%s_fix_clado.txt' % strain, 'r')
for line in f:
	
	sample_bam_name = line.strip().split('.ready.bam')[0].split('.bam')[0]
	sample_name = bam_name_to_clone_name(sample_bam_name)
	
	# Compound
	if sample_bam_name not in clone_compound_dict:
		print("Skipping %s" % sample_name)
		continue
	
	if clone_strain_dict[sample_bam_name] != strain: # Wrong strain
		print("WRONG STRAIN FOR %s" % sample_name)
		continue
	
	compound = clone_compound_dict[sample_bam_name]
	
	# Get coverage info
	total_cov = sample_coverage_dict[sample_name]
	avg_cov = sample_avg_coverage_dict[sample_name]
	perc_bases_gt_20 = sample_perc_bases_gt_20_dict[sample_name]
	
	# Get copy ratio info
	copyratio_data_available = True
	try:
		chrom_interval_copyratio_dict = get_copyratio_dict(sample_bam_name, strain)
	except:
		copyratio_data_available = False
	
	# Get predicted CNVs
	predicted_CNVs = get_predicted_CNVs(sample_bam_name)
	if len(predicted_CNVs) == 0:
		print("No CNVs found for %s" % sample_name)
		continue
	
	# Combine adjacent CNVs
	CNV_chrom_start_end_dict = defaultdict(dict); CNV_chrom_end_start_dict = defaultdict(dict)
	for chrom, start, end in predicted_CNVs:
		CNV_chrom_start_end_dict[chrom][start] = end; CNV_chrom_end_start_dict[chrom][end] = start
	
	for chrom in CNV_chrom_start_end_dict:
		for start1 in CNV_chrom_start_end_dict[chrom].keys():
			end1 = CNV_chrom_start_end_dict[chrom][start1]
			for end2 in CNV_chrom_end_start_dict[chrom].keys():
				start2 = CNV_chrom_end_start_dict[chrom][end2]
				if (start1 - end2) > 0 and (start1 - end2) < 40000:
					predicted_CNVs.append((chrom, start2, end1))
	
	# Get reads that may serve as evidence for amplification
	chrom_reads_dict = get_supporting_reads_dict(sample_bam_name)
	if not chrom_reads_dict:
		print("No supporting read file for %s" % sample_bam_name)
		continue
	
	# Get coverage information
	num_RF_chrom = sum([len(chrom_reads_dict['RF'][chrom]) for chrom in chrom_reads_dict['RF']])
	num_RR_chrom = sum([len(chrom_reads_dict['RR'][chrom]) for chrom in chrom_reads_dict['RR']])
	num_FF_chrom = sum([len(chrom_reads_dict['FF'][chrom]) for chrom in chrom_reads_dict['FF']])
	
	avg_RF_prop_chrom = num_RF_chrom/float(total_cov) if total_cov != 0 else -1
	avg_RRFF_prop_chrom = (num_RR_chrom + num_FF_chrom)/float(total_cov) if total_cov != 0 else -1
	
	# print("Working on %s..." % sample_name)
	
	for CNV in predicted_CNVs:
		CNV_chrom, CNV_start, CNV_end = CNV
		
		# Consider deletion
		if copyratio_data_available:
			copyratio_list = get_copyratio_list(chrom_interval_copyratio_dict, CNV_chrom, CNV_start, CNV_end)
			average_copyratio = np.mean(copyratio_list)
		else:
			average_copyratio = 0 # Unknown
		
		if average_copyratio < 0:
			classification = "POTENTIAL DELETION"; new_start = None; new_end = None; support = None
		else:
			classification, new_start, new_end, support = validate_CNV(CNV, chrom_reads_dict)
			if type(support) == 'tuple':
				support = min(support)
		
		if (sample_name, CNV_chrom, CNV_start, CNV_end) in CNV_pval_dict:
			pval = CNV_pval_dict[(sample_name, CNV_chrom, CNV_start, CNV_end)]
		else:
			pval = 'combined'
		
		items = [compound, sample_name, CNV_chrom, CNV_start, CNV_end, pval, average_copyratio, classification, support, new_start, new_end, avg_cov, perc_bases_gt_20, avg_RF_prop_chrom, avg_RRFF_prop_chrom]
		output_f.write('\t'.join([str(item) for item in items]) + '\n')

f.close()
output_f.close()	
	
