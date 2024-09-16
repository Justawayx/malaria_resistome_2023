from collections import defaultdict
import pickle
import sys

# =================
# Parameters
# =================

index = int(sys.argv[1])
chromosomes = ['Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3', 'Pf3D7_04_v3', 'Pf3D7_05_v3', 'Pf3D7_06_v3', 'Pf3D7_07_v3', 'Pf3D7_08_v3', 'Pf3D7_09_v3', 'Pf3D7_10_v3', 'Pf3D7_11_v3', 'Pf3D7_12_v3', 'Pf3D7_13_v3', 'Pf3D7_14_v3']
chrom = chromosomes[index]

output_dir = '/projects/winzeler/ROTATION_PROJECT/daisy/dnds/pickles'

alt_depth_threshold = 50
alt_freq_threshold = 0.5

# =================
# Helper functions
# =================

def convert_format_dict(format_str, format_items_str):
    format_cats = format_str.split(':')
    format_items = format_items_str.split(':')
    return {cat: item for cat, item in zip(format_cats, format_items)}

def convert_info_dict(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        parts = item.split('=')
        if len(parts) == 2:
            cat, val = parts
        else:
            cat = parts[0]; val = ''
        info_dict[cat] = val
    return info_dict

# =================
# Store Pf6 SNVs
# =================

print("Working on %s" % chrom)

f = open('/projects/winzeler/ROTATION_PROJECT/daisy/pf6/Pf_60_public_%s.ann.txt' % chrom, 'r')

# Skip the meta-information lines
for line in f:
    if line[:2] !=  '##':
        break

samples = line.strip().split('\t')[9:] # Header

sample_chrom_pos_syn_snv_dict = {s: defaultdict(dict) for s in samples} # sample -> chrom -> pos -> alt allele
sample_chrom_pos_nonsyn_snv_dict = {s: defaultdict(dict) for s in samples} # sample -> chrom -> pos -> alt allele
chrom_gene_pos_dict = {c: defaultdict(set) for c in chromosomes} # chrom -> gene ID -> pos with variants in the above dicts

for line in f:
    items = line.strip().split('\t')
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
    
    if FILTER != 'PASS' or QUAL == '.' or float(QUAL) < 500 or (REF not in ['A', 'G', 'C', 'T']) or ('SNPEFF_TRANSCRIPT_ID' not in INFO): # Filter
        continue
    
    pos = int(POS)    
    alt_alleles = ALT.split(',')
    info_items = convert_info_dict(INFO)
    variant_types = [eff_str.split('(')[0] for eff_str in info_items['EFF'].split(',')]
    gene_id = info_items['SNPEFF_TRANSCRIPT_ID']
    
    for tup, sample in zip(items[9:], samples):
        FITEMS = convert_format_dict(FORMAT, tup)
        if FITEMS['GT'] == './.' or FITEMS['DP'] == '.': # Filter per sample
            continue
        
        alt_allele_depths = [int(val) for val in FITEMS['AD'].split(',')[1:]]; total_depth = int(FITEMS['DP'])
        
        for alt_allele, alt_depth, variant_type in zip(alt_alleles, alt_allele_depths, variant_types):
            if alt_allele in ['A', 'G', 'C', 'T'] and alt_depth >= alt_depth_threshold and (alt_depth/total_depth) >= alt_freq_threshold: # Confident SNV
                
                if 'synonymous_variant' in variant_type:
                    sample_chrom_pos_syn_snv_dict[sample][CHROM][pos] = alt_allele
                    chrom_gene_pos_dict[CHROM][gene_id].add(pos)
                elif 'missense_variant' in variant_type:
                    sample_chrom_pos_nonsyn_snv_dict[sample][CHROM][pos] = alt_allele
                    chrom_gene_pos_dict[CHROM][gene_id].add(pos)    

print("Done!")

# =================
# Pickle
# =================

pickle.dump(sample_chrom_pos_syn_snv_dict, open('%s/pf6_%s_sample_chrom_pos_syn_snv_dict.pkl' % (output_dir, chrom), 'wb'))
pickle.dump(sample_chrom_pos_nonsyn_snv_dict, open('%s/pf6_%s_sample_chrom_pos_nonsyn_snv_dict.pkl' % (output_dir, chrom), 'wb'))
pickle.dump(chrom_gene_pos_dict, open('%s/pf6_%s_chrom_gene_pos_dict.pkl' % (output_dir, chrom), 'wb'))
