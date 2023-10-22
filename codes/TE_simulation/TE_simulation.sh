#!/bin/bash

# This bash script aims to generate pseudo numts by reshuffling their genomic coordinates while maintaining their length and estimate the percentage of TE in their flanking regions (up to 5kb) in a 500bp window. This process is repeated 1000 times and will be used as background in comparison with the observed data (the percentage of TE in the flanking regions of real numts.

# This script uses 'cattle' as an example. To run the script, you need:

# "Cattle_MB_numts.bed": NUMT genomic coordinates of the cattle;
# "Btau_length": The length of scaffolds in the Bos taurus genome;
# "Btau_excl": The genomic coordinates in bed format that are excluded from reshuffling numts coordinates at the end of scaffolds (5kb at the beginning and 5kb at the end); this is to avoid reshuffling numts close to the end of scaffolds so that the flanking regions would be less than 5kb.
# "length.pl": a perl script to calculate the length of scaffolds in a genome fasta file.


# Repeat the code for 1000 times;
for ((i = 1; i <= 1000; i++)); do

# Reshuffle cattle numt genomic coordinates;
bedtools shuffle -i Cattle_MB_numts.bed -g Btau_length -excl Btau_excl >shuffled.bed

# Obtain the 5kb upstream and downstream flanking regions in bed format;
bedtools flank -i shuffled.bed -g Btau_length -b 5000 >flank.bed

# Separate the flanking regions to 5kb upstream and downstream, respectively;
awk 'NR % 2 == 1' flank.bed > Btau_upstream.bed
awk 'NR % 2 == 0' flank.bed > Btau_downstream.bed

# Obtain the 5kb sequences for upstream and downstream flanking regions;
bedtools getfasta -bed Btau_upstream.bed -fi GCF_002263795.1_ARS-UCD1.2_genomic.fna.masked -fo Btau_upstream.fa
bedtools getfasta -bed Btau_downstream.bed -fi GCF_002263795.1_ARS-UCD1.2_genomic.fna.masked -fo Btau_downstream.fa

# Obtain the length of each sequence; (This is required to obtain sequences in windows.)
perl length.pl Btau_downstream.fa test_down
perl length.pl Btau_upstream.fa test_up

# Chop the 5kb upstream and downstream flanking regions into 500bp windows (10 windows each)
bedtools makewindows -g test_up -w 500 > up_windows.bed
bedtools getfasta -fi Btau_upstream.fa -bed up_windows.bed -fo up_window.fa
bedtools makewindows -g test_down -w 500 > down_windows.bed
bedtools getfasta -fi Btau_downstream.fa -bed down_windows.bed -fo down_window.fa

# Estimate average TE percentage per window by calculating masked regions (N) for upstream and downstream flanking regions respectively; 
awk '/^>/ {if (seq) print count_N / length(seq) * 100; seq=""; count_N=0; print; next} {seq = seq $0} /N/ {count_N += gsub(/N/, "")} END {if (seq) print count_N / length(seq) * 100;}' down_window.fa >down_perc
awk '/^>/ {if (seq) print count_N / length(seq) * 100; seq=""; count_N=0; print; next} {seq = seq $0} /N/ {count_N += gsub(/N/, "")} END {if (seq) print count_N / length(seq) * 100;}' up_window.fa >up_perc

# Write the average TE percentage of upstream flanking regions per window into a file;
awk 'NR >= 2 && NR % 20 == 2' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 4 && NR % 20 == 4' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 6 && NR % 20 == 6' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 8 && NR % 20 == 8' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 10 && NR % 20 == 10' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 12 && NR % 20 == 12' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 14 && NR % 20 == 14' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 16 && NR % 20 == 16' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 18 && NR % 20 == 18' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_up.txt
awk 'NR >= 20 && NR % 20 == 0' up_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\n", sum / count; }' >> result_up.txt

# Write the average TE percentage of downstream flanking regions per window into a file;
awk 'NR >= 2 && NR % 20 == 2' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 4 && NR % 20 == 4' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 6 && NR % 20 == 6' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 8 && NR % 20 == 8' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 10 && NR % 20 == 10' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 12 && NR % 20 == 12' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 14 && NR % 20 == 14' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 16 && NR % 20 == 16' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 18 && NR % 20 == 18' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\t", sum / count; }' >> result_down.txt
awk 'NR >= 20 && NR % 20 == 0' down_perc | awk '{ sum += $1; count++ } END { if (count > 0) printf "%f\n", sum / count; }' >> result_down.txt

# Remove intermediate files. The two files (result_up.txt; result_down.txt) contains the average TE percentage per window (10 windows in total) respectively;
rm Btau_downstream.* Btau_upstream.* down* shuffled.bed test* up*

  echo "Iteration $i"
done
