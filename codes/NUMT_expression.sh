#!/bin/bash

# This bash script aims to evaluate expression levels of NUMTs. To achieve this, we obtained each NUMT sequence along with their 200bp upstream and downstream flanking regions for each species, and quantified their expression using publicly available RNA-Seq datasets. The NUMTs that are deemed to be expressed should satisfy the following criteria: 1) at least 2 RNA-Seq reads support the junctions between NUMTs and their flanking regions with at least 5bp overhangs; 2) the coverage of NUMTs by RNA-Seq reads is over 70%.

# This script uses the human brain sample (SRR8942877) as an example. You will need the following files:
# Raw RNA-Seq data: SRR8942877.fastq
# The human genome file: GCF_000001405.39_GRCh38.p13_genomic.fna
# The length of each scaffold in the human genome file: Hsap_length.txt
# The genomic coordinates of human NUMTs in .bed format: Human_MB_numts.bed
# The length of each human NUMT in .bed format: Numt_length.bed


# Obtain human NUMT sequences as well as their 200bp upstream and downstream flanking regions:

bedtools slop -i Human_MB_numts.bed -g Hsap_length.txt -b 200 >Hsap_numt_ref.fa


# Quality control of raw RNA-Seq data:

cutadapt -m 75 -q 25,25 -o HQ_SRR8942877.fastq SRR8942877.fastq


# Index reference sequences:

hisat2-build Hsap_numt_ref.fa ref


# Map preprocessed RNA-Seq reads to the references:

hisat2 -x ref -U HQ_SRR8942877.fastq -S Mapped.sam


# Convert and sort the SAM file to the Bed file:

samtools view -S -b Mapped.sam >Mapped.bam
samtools sort Mapped.bam >Mapped.sorted.bam
samtools index Mapped.sorted.bam
bedtools bamtobed -i Mapped.sorted.bam Mapped.sorted.bam >Mapped.bed


# Estimate the coverage of RNA-Seq reads for each NUMT:

bedtools coverage -a Numt_length.bed -b Mapped.bed >Mapped_coverage.txt


# The Mapped.bam file could be visualised using the IGV genome browser.