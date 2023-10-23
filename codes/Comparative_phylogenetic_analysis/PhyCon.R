# This script is to estimate phylogenetic correlation between two traits.
# The time-calibrated phylogenetic tree and the trait table are provided in the same folder as this R script.

library(phytools)

# the time-calibrated phylogenetic tree of 45 mammals
tree=read.tree("Timetree_all.tre")

# the file that contains genome statistics (e.g. genome size (gz), the number of high score pairs (hsp))
table=read.table("genome_stats.txt",header=TRUE)

# extract two variables for the analysis: 'gz' and 'the number of hsp' as examples
X=cbind(table$gz,table$hsp)

# rename the variables
colnames(X)=c("gz","hsp")

# establish the covariance matrix
obj<-phyl.vcv(X,vcv(tree),1)

# phylogenetic correlation between 'gz' and 'hsp'
r.xy=cov2cor(obj$R)["gz","hsp"]
r.xy

# t-statistics & P-value
t.xy<-r.xy*sqrt((Ntip(tree)-2)/(1-r.xy^2))
P.xy<-2*pt(abs(t.xy),df=Ntip(tree)-2,lower.tail=F)
P.xy
