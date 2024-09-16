import copy
import fractions

nucs = ['A', 'G', 'C', 'T']

alt_nuc_dict = {nuc: [altnuc for altnuc in nucs if (altnuc != nuc)] for nuc in nucs}

f = open("/projects/winzeler/ROTATION_PROJECT/daisy/genomics_utils/utils/codons.txt", 'r')
header = f.readline()

codon_aa_dict = {}
codon_aa1_dict = {}

for line in f:
	codon, aa_3, aa_1, aa_F = line.strip().split('\t')
	codon_aa_dict[codon] = aa_3
	codon_aa1_dict[codon] = aa_1

def classify_snv(nuc_sequence, pos, new_allele, verbose=False, formatted=False):
	codon_start_pos = 3*(pos // 3)
	relative_pos = pos % 3
	codon = nuc_sequence[codon_start_pos:codon_start_pos+3]
	old_aa = codon_aa_dict[codon]
	
	lcodon = list(codon)
	new_lcodon = copy.copy(lcodon); new_lcodon[relative_pos] = new_allele
	new_codon = ''.join(new_lcodon)
	new_aa = codon_aa_dict[new_codon]
	
	if new_codon == codon:
		classification = 'no_change'
	elif new_aa == old_aa:
		classification = "synonymous"
	else:
		if new_aa == 'Stp':
			classification = "nonsense"
		else:
			classification = "nonsynonymous"
	
	# Example: p.Leu322Phe/c.964C>T
	if verbose and formatted:
		return (classification, ("%i%s>%s" % (pos+1, lcodon[relative_pos], new_allele)), ("%s%i%s" % (old_aa, (codon_start_pos//3) + 1, new_aa)))
	elif verbose:
		return (classification, (codon, new_codon), (old_aa, new_aa))
	else:
		return classification

def frac_nonsyn(codon, pos):
	nuc = codon[pos]
	alt_nucs = alt_nuc_dict[nuc]
	lcodon = list(codon)
	old_aa = codon_aa_dict[codon]
	
	num_nonsyn = 0
	for alt_nuc in alt_nucs:
		new_lcodon = copy.copy(lcodon); new_lcodon[pos] = alt_nuc
		new_aa = codon_aa_dict[''.join(new_lcodon)]
		if new_aa != old_aa:
			num_nonsyn += 1
	
	return fractions.Fraction(num_nonsyn, len(alt_nucs))

def codon_to_expected_n_s(codon):
	n = 0
	for pos in range(3):
		n += frac_nonsyn(codon, pos)
	s = 3 - n
	return (n, s)

def nuc_sequence_to_expected_n_s(nuc_sequence):
	total_n = 0; total_s = 0
	for i in range(0, len(nuc_sequence), 3):
		codon = nuc_sequence[i:i+3]
		n, s = codon_to_expected_n_s(codon)
		total_n += n
		total_s += s
	return (float(total_n), float(total_s))

def nuc_sequence_to_aa_sequence(nuc_sequence):
	aa_sequence = ''
	for i in range(0, len(nuc_sequence), 3):
		codon = nuc_sequence[i:i+3]
		aa = codon_aa1_dict[codon]
		if aa == 'O':
			break
		aa_sequence += aa
	return aa_sequence

# Testing
nuc_sequence = '''ATGATGAGAAAATTAGCTATTTTATCTGTTTCTTCCTTTTTATTTGTTGAGGCCTTATTCCAGGAATACCAGTGCTATGGAAGTTCGTCAAACACAAGGGTTCTAAATGAATTAAATTATGATAATGCAGGCACTAATTTATATAATGAATTAGAAATGAATTATTATGGGAAACAGGAAAATTGGTATAGTCTTAAAAAAAATAGTAGATCACTTGGAGAAAATGATGATGGAAATAACGAAGACAACGAGAAATTAAGGAAACCAAAACATAAAAAATTAAAGCAACCAGCGGATGGTAATCCTGATCCAAATGCAAACCCAAATGTAGATCCCAATGCCAACCCAAATGTAGATCCAAATGCAAACCCAAATGTAGATCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCCAATGCAAATCCTAATGCAAACCCAAATGCAAACCCAAACGTAGATCCTAATGCAAATCCAAATGCAAACCCAAACGCAAACCCCAATGCAAATCCTAATGCAAACCCCAATGCAAATCCTAATGCAAATCCTAATGCCAATCCAAATGCAAATCCAAATGCAAACCCAAACGCAAACCCCAATGCAAATCCTAATGCCAATCCAAATGCAAATCCAAATGCAAACCCAAATGCAAACCCAAATGCAAACCCCAATGCAAATCCTAATAAAAACAATCAAGGTAATGGACAAGGTCACAATATGCCAAATGACCCAAACCGAAATGTAGATGAAAATGCTAATGCCAACAGTGCTGTAAAAAATAATAATAACGAAGAACCAAGTGATAAGCACATAAAAGAATATTTAAACAAAATACAAAATTCTCTTTCAACTGAATGGTCCCCATGTAGTGTAACTTGTGGAAATGGTATTCAAGTTAGAATAAAGCCTGGCTCTGCTAATAAACCTAAAGACGAATTAGATTATGCAAATGATATTGAAAAAAAAATTTGTAAAATGGAAAAATGTTCCAGTGTGTTTAATGTCGTAAATAGTTCAATAGGATTAATAATGGTATTATCCTTCTTGTTCCTTAATTAG'''

test_aa_sequence = nuc_sequence_to_aa_sequence(nuc_sequence)

aa_sequence = '''MMRKLAILSVSSFLFVEALFQEYQCYGSSSNTRVLNELNYDNAGTNLYNELEMNYYGKQENWYSLKKNSRSLGENDDGNNEDNEKLRKPKHKKLKQPADGNPDPNANPNVDPNANPNVDPNANPNVDPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNVDPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNKNNQGNGQGHNMPNDPNRNVDENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLIMVLSFLFLN'''

print(aa_sequence == test_aa_sequence)
exp_n, exp_s = nuc_sequence_to_expected_n_s(nuc_sequence)
print(exp_n, exp_s)
print(exp_n + exp_s, len(nuc_sequence))
