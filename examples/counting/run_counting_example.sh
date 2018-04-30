# Run counting on fastqs

# define samples
samples="A549-CV4-1000_d21_1 A549-CV4-1000_d21_2 A549-CV4-1000_d28_1 A549-CV4-1000_d28_2"

# make output directories
mkdir -p counts
mkdir -p alignments

for sample in $samples;
do 
	ctg count --config count_config.txt \
		--sample ${sample} \
		-1 fastqs/${sample}_R1.fastq \
		-2 fastqs/${sample}_R2.fastq \
		--output_counts counts/${sample}.counts.csv \
		--output_bam alignments/${sample}.alignments.bam
done
