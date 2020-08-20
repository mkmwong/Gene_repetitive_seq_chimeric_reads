#!/bin/bash
#SBATCH --job-name=hisat2_chime_reads
#SBATCH --time=02:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH -p hns,normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkmwong@stanford.edu

###################################################################################
## !! This should be the only session to make changes to !! ##
## Software requirement: hisat2, samtools, bowtie2, R 
## Package requirement in R: GenomicRanges, biomaRt, dplyr
## Currently only support paired-end RNA-Seq file
## Genome index should be created with hisat2
## Repeat index should be created with bowtie2
input_file_path="/scratch/groups/ashbym/mandy/novogene/raw_data/RbKO2"
fastq_1="RbKO2_FRRL202355629-1a_H3TG5DSXY_L3_1.fq.gz"
fastq_2="RbKO2_FRRL202355629-1a_H3TG5DSXY_L3_2.fq.gz"
out_path="/scratch/groups/ashbym/mandy/novogene/raw_data/RbKO1/pipeline_out"
hisat2_index_path="/home/groups/ashbym/mandy/hisat2_index/hg19_hisat2"
out_file_name="RbKO1_hisat.sam"
install_dir="/scratch/groups/ashbym/mandy/novogene/raw_data"
repeat_genome=$install_dir"/repeatome/active_rep"
###################################################################################

mkdir $out_path
cd $out_path
module load R

## Align fastq file with hisat2
echo "Aligning reads using hisat2......"
hisat2 -x $hisat2_index_path -1 $input_file_path"/"$fastq_1 -2 $input_file_path"/"$fastq_2 -S $out_path"/"$out_file_name

## Find reads with mate unmapped and unmapped reads with mate mapped
echo "Seperating and filtering reads for aligning to repeatome......"
samtools view -bS -h -f 8 -F 4  $out_file_name > mapped_reads.bam
samtools view -bS -h  -F 8 -f 4 $out_file_name > unmapped_mate.bam

## Removing secondary alignment of mapped reads
samtools view -bh -F 256 mapped_reads.bam > mapped_reads_nosec.bam

## Separating reads based on their orientation and R1/R2
samtools view -b -h -F 16 -f 64 mapped_reads_nosec.bam > mapped_reads_forward_R1.bam
samtools view -b -h -F 16 -f 128 mapped_reads_nosec.bam > mapped_reads_forward_R2.bam
samtools view -b -h -f 16 -f 64 mapped_reads_nosec.bam > mapped_reads_reverse_R1.bam
samtools view -b -h -f 16 -f 128 mapped_reads_nosec.bam > mapped_reads_reverse_R2.bam

## Convert unmapped reads back to fastq
samtools fastq -1 unmapped_1.fastq -2 unmapped_2.fastq -s unmapped_singlet.fastq unmapped_mate.bam

## Mapping unmapped reads to "repeatome"
echo "Aligning reads to repeatome......"
bowtie2 -x $repeat_genome unmapped_singlet.fastq | samtools view -bS > repeat_out.bam

## Only retain unmapped reads that are mapped to repeat
samtools view -F 4 repeat_out.bam > mapped_to_repeat1.sam

# Isolate readnames of reads from step above
cut -f1 mapped_to_repeat1.sam > readname

## Filter isolated mapped reads by read names
samtools view  mapped_reads_forward_R1.bam | fgrep -w -f readname > transcript_forward_R1.sam
samtools view  mapped_reads_reverse_R1.bam | fgrep -w -f readname > transcript_reverse_R1.sam
samtools view  mapped_reads_forward_R2.bam | fgrep -w -f readname > transcript_forward_R2.sam
samtools view  mapped_reads_reverse_R2.bam | fgrep -w -f readname > transcript_reverse_R2.sam

## Subset sam files to retain only the useful columns
cut -f1,3,4,6 transcript_forward_R1.sam | sort -k 1 > gene_sub_for_R1
cut -f1,3,4,6 transcript_reverse_R1.sam | sort -k 1 > gene_sub_rev_R1
cut -f1,3,4,6 transcript_forward_R2.sam | sort -k 1 > gene_sub_for_R2
cut -f1,3,4,6 transcript_reverse_R2.sam | sort -k 1 > gene_sub_rev_R2
cut -f1,3,4,6 mapped_to_repeat1.sam | sort -k 1 > repeat_sub

## Join the mapped with with the unmapped reads mapped to repeat
join gene_sub_for_R1 repeat_sub > coor_for_R1.out
join gene_sub_rev_R1 repeat_sub > coor_rev_R1.out
join gene_sub_for_R2 repeat_sub > coor_for_R2.out
join gene_sub_rev_R2 repeat_sub > coor_rev_R2.out

## run RScript for removing the reads that are just natually existing repeats
Rscript $install_dir"/filter_reads.R" $out_path $install_dir

