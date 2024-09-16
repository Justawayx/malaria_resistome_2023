import math, random
from collections import defaultdict
import pickle, gzip, csv
from urllib.parse import unquote

from utils.dnds_utils import (classify_snv, nuc_sequence_to_expected_n_s, nuc_sequence_to_aa_sequence)

from Bio import SeqIO

comp_dict = {'G': 'C', 'C': 'G', 'T': 'A', 'A': 'T', 'g': 'c', 'c': 'g', 't': 'a', 'a': 't'}

def reverse_complement(dna):
	complement = ''
	for nuc in dna[::-1]:
		complement += comp_dict[nuc]
	return complement

# ============================
# Load genome
# ============================

fasta_fpath = "/projects/winzeler/GENOME_RESOURCES/pf/p_fal_ref/p_fal.fasta"

chromosome_sequence_dict = {}

for seq_record in SeqIO.parse(fasta_fpath, "fasta"):
	chrom = str(seq_record.id)
	sequence = str(seq_record.seq)
	chromosome_sequence_dict[chrom] = sequence

# ============================
# Load gene information
# ============================

gff_fpath = "/projects/winzeler/GENOME_RESOURCES/pf/p_fal_ref/p_fal.gff"

chromosomes = []
for line in open(gff_fpath, 'r'):
	if line[0] != '#':
		break
	if line[:17] == '##sequence-region':
		chromosomes.append(line.strip().split()[1])

chrom_gene_exon_interval_dict = {chrom: defaultdict(dict) for chrom in chromosomes} # chrom -> gene ID -> exon ID -> (start, end, strand_direction)

chrom_gene_ids_dict = {chrom: set() for chrom in chromosomes} # All protein coding gene IDs

chrom_gene_desc_dict = {chrom: {} for chrom in chromosomes} # gene_id -> description
chrom_gene_interval_dict = {chrom: {} for chrom in chromosomes} # gene_id -> (start, end, strand_direction)

for line in open(gff_fpath, 'r'):
	if line.strip() == '##FASTA':
		break
	if line[0] == '#':
		continue
	chrom, source, feature_type, start_pos, end_pos, _, sdirection, _, info = line.strip().split('\t')
	start_pos = int(start_pos); end_pos = int(end_pos)
	if feature_type == 'exon':
		exon_id = info.split(';')[0].split('ID=')[1]
		gene_id = exon_id.split('exon_')[1].split('-')[0]
		chrom_gene_exon_interval_dict[chrom][gene_id][exon_id] = (start_pos, end_pos, sdirection)
	if feature_type == 'gene':
		gene_id = info.split(';')[0].split('ID=')[1]
		gene_desc = info.split(';')[2].split('description=')[1]
		chrom_gene_desc_dict[chrom][gene_id] = unquote(gene_desc.strip()).replace('+', ' ')
		chrom_gene_interval_dict[chrom][gene_id] = (start_pos, end_pos, sdirection)
	if feature_type == 'CDS': # Actually coding
		gene_id = info.split(';')[0].split('ID=')[1].split('cds_')[1].split('-')[0]
		chrom_gene_ids_dict[chrom].add(gene_id)

# ================================================
# Compute expected N and expected S for each gene
# ================================================

gene_expn_exps_dict = {}

for chrom in chromosomes:
	print("Working on %s..." % chrom)
	sequence = chromosome_sequence_dict[chrom]
	
	for gene_id in chrom_gene_ids_dict[chrom]:
		
		# Weird exception
		if gene_id == 'PF3D7_0112400':
			continue
		
		# Ignore pseudogenes
		if 'pseudogene' in chrom_gene_desc_dict[chrom][gene_id]:
			continue
		
		forward_nuc_sequence = ''
		
		exon_ids = list(chrom_gene_exon_interval_dict[chrom][gene_id].keys())
		sorted_exon_ids = sorted(exon_ids, key=lambda x: int(x.split('-')[1]))
		for exon_id in sorted_exon_ids:
			start, end, sdirection = chrom_gene_exon_interval_dict[chrom][gene_id][exon_id]
			genome_nuc_sequence = sequence[start-1:end]
			if sdirection == '-':
				exon_forward_nuc_sequence = reverse_complement(genome_nuc_sequence)
				forward_nuc_sequence += exon_forward_nuc_sequence
			elif sdirection == '+':
				exon_forward_nuc_sequence = genome_nuc_sequence
				forward_nuc_sequence += exon_forward_nuc_sequence
		
		expn, exps = nuc_sequence_to_expected_n_s(forward_nuc_sequence)
		gene_expn_exps_dict[gene_id] = (expn, exps)

# ================================================
# Store data in table
# ================================================

f = open("/projects/winzeler/ROTATION_PROJECT/daisy/pfal_expected_ns_s.csv", 'w')
f.write(','.join(['Chromosome', 'Gene ID', 'Expected NS', 'Expected S']) + '\n')

for chrom in chromosomes[2:] + chromosomes[:2]:
	for gene_id in chrom_gene_ids_dict[chrom]:
		
		# Weird exception
		if gene_id == 'PF3D7_0112400':
			continue
		
		# Ignore pseudogenes
		if 'pseudogene' in chrom_gene_desc_dict[chrom][gene_id]:
			continue
		
		expn, exps = gene_expn_exps_dict[gene_id]
		f.write(','.join([str(val) for val in [chrom, gene_id, expn, exps]]) + '\n')

f.close()
