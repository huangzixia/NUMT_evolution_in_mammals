# The script "genenral_commands.sh"

## 1) Assessment of genome quality using BUSCO
## 2) Identification of NUMTs from genomes using BLAST
## 3) Masking transposable elements (TE) using repeatmasker
## 4) Authentication of NUMTs with PacBio raw reads using BLAST
## 5) Obtaining ten up- and down-stream protein-coding genes to infer NUMT genomic synteny

# The folder "Comparative_phylogenetic_analysis"

## 1) The R code (PhyCon.R) evaluates the phylogenetic correlation between two variables
## 2) The phylogenetic tree (Timetree_all.tre) includes the time-calibrated phylogenetic tree across 45 mammals. This is needed for the phylogenetic correlation analysis
## 3) The file (genome_stats.txt) listed all the genomic traits investigated in this study

# The folder "NUMT_expression"

## 1) The file (NUMT_expression.sh) documented the commands used to perform the expression analysis of NUMTs
## 2) The file (Human_MB_numts.bed) includes the genomic coordinates of NUMTs/NUMT blocks in the human genome
## 3) The file (Hsap_numt_ref.fa) includes the genomic sequences of NUMTs/NUMT blocks as well as their 200 bp flanking regions in both ends
## 4) The file (Numt_length.bed) contains the mtDNA coordinates of NUMTs/NUMT blocks. This is required for evaluating the RNA-Seq read coverage for each NUMT
## 5) The file (Hsap_length.txt) contains the information of scaffold length in the human genome

# The folder "Coverage_simulation" 

## 1) The R code (Coverage_simulation.R) aims to generate the simulation data of NUMT mtDNA coverage and perform pairwise comparisons across windows with a size of 50 bp
## 2) The file (NUMT_length_species.txt) contains the length information of each NUMT per species, which is required for simulating data

# The folder "TE_simulation"

## 1) The script (TE_simmulation.sh) aims to generate pseudo-numts and estimate the percentage of TE in their flanking regions which will serve as background data for comparisons with that of the observed data
## 2) The script (length.pl) is to obtain the sequence length which is required to run the TE_simulation.sh script
## 3) The file (Cattle_MB_numts.bed) contains the genomic coordinates of NUMTs/NUMT blocks in the cattle genome
## 4) The file (Btau_length) contains the length of each scaffold in the cattle genome
## 5) The file (Btau_excl) contains the genomic coordinates that are excluded when reshuffling NUMT genomic coordinates.
