bam_dir=/storage/NFS/ANALYSIS/DNAseq/Madeline/BAMS/missing_2023/3D7
bam_filelist=/storage/NFS/ANALYSIS/DNAseq/test/bamnames_3D7_fix.txt
ref_fasta=/storage/NFS/GENOME_RESOURCES/pf/p_fal_ref/p_fal.fasta
ref_map=/storage/NFS/GENOME_RESOURCES/pf/p_fal_ref/map.fa.gz
output_dir=/storage/NFS/ANALYSIS/DNAseq/test/delly_output_3D7_fix

for bamname in $(cat $bam_filelist)
do
	x=${bamname%.ready.bam}
	y=${x##*/}
	z=${y%.bam}
	sample=${z##*/}
	echo $sample
	delly cnv -g $ref_fasta -m $ref_map $bam_dir/$bamname > $output_dir/${sample}.vcf
done

