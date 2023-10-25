# This file includes the commands that were used to mask genomes, assess genome quality, identify NUMTs from genomes, as well as
# other general data analyses. The scripts for specific data analyses, e.g. NUMT coverage simulation, comparative phylogenetic
# correlation analysis, could be found in separate folders in Github at:
# https://github.com/huangzixia/NUMT_evolution_in_mammals/tree/main/codes
# Please replace the files in brackets [] with your own files.


# 1. Assessing genome quality using BUSCO (v4.0)

busco -m genome -i [genome_file.fa] -l eukaryota_odb10 -o [output_name] 


# 2.Identify NUMTs from genomes using dc-megablast with optimal parameters (BLAST v2.9.0). Other BLAST methods (e.g. megablast and blastn) and parameters used to optimise the pipeline can be found in Supplementary Table S10

## Make BLAST database

makeblastdb -in [genome_file.fa] -dbtype nucl -out [genome_index_name]

## Run dc-megablast with the optimal parameters

blastn -db [genome_index_name] -query [mtDNA_sequence.fa] -out [blast_output] -evalue e-3 -template_type optimal -template_length 18 -word_size 11 -outfmt 7 -name


# 3. Mask genomes using repeatmasker (v4.1.2) for transposable element (TE) analyses of NUMT flanking regions. 

## Obtain the 5kb upstream and downstream flanking regions of NUMTs per species

bedtools flank -i [NUMT_genomic_coordinates.bed] -g [scaffold_length.txt] -b 5000 >[flank.bed]

## Obtain 5kb flanking sequences of NUMTs

bedtools getfasta -fi [genome_file.fa] -bed [flank.bed] -fo [flank.fa]

## Run repeatmasker 

RepeatMasker -pa 24 -norna -dir [directory] [flank.fa]


# 4. Authentication of predicted NUMTs using optimal dc-megablast parameters aforementioned, with NUMT sequences as query and PacBio raw reads as database

## Make BLAST database

makeblastdb -in [PacBio_raw_reads.fa] -dbtype nucl -out [genome_index_name]

## Run dc-megablast with the optimal parameters

blastn -db [genome_index_name] -query [mtDNA_sequence.fa] -out [blast_output] -evalue e-3 -template_type optimal -template_length 18 -word_size 11 -outfmt 7 -name


# 5. Obtain ten up- and down-stream protein-coding genes (anchors) of each NUMT per species for synteny analyses

## Convert GFF annotation files to BED files

convert2bed -i gff < [genome_annotation.gff] > [genome_annotation.bed]

## Filter the genome annotation files. We only retained the feature of protein-coding genes in the annotation files

sed '/exon/d' [genome_annotation.bed] > [filtered_exon.bed]
sed '/CDS/d' [filtered_exon.bed] > [filtered_exon_cds.bed]
sed '/rna/d' [filtered_exon_cds.bed] > [filtered_rna.bed]
sed '/LOC/d' [filtered_rna.bed] > [filtered_rna2.bed]
sed '/x[0-9]\{2,\}/d' [filtered_rna2.bed] > [filtered_rna2_no.bed]
sed '/x[2-9]/d' [filtered_rna2_no.bed] > [filtered_noIso.bed]
sed '/region/d' [filtered_noIso.bed] > [filtered2_noIso2.bed]
sed 's/gene-//g' -i [filtered2_noIso2.bed]
sed 's/match//g' -i [filtered2_noIso2.bed]
bedtools sort -i [filtered2_noIso2.bed] > [filtered2_noIso2.sort.bed]
bedtools merge -i [filtered2_noIso2.sort.bed] -c 4 -o distinct > [filtered2_noIso2.merged.bed]

## Merge genomic coordinates of NUMTs with protein-coding genomic coordinates

cat [Hsap_MB_NUMT.bed] [filtered2_noIso2.merged.bed] > [merged.bed]
bedtools sort -i [merged.bed] > [merged.sort.bed]
bedtools merge -i [merged.sort.bed] -c -o distinct > [synteny.bed]

## Obtain 10 protein-coding genes upstream and 10 genes downstream of each NUMT per species

grep -C 10 numt [syteny.bed] >synteny_list.txt

## Change the file format

pr [syteny_list.txt] -22 -a -t -s >[Hsap_syteny.txt]





