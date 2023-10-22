# This R script is to generate the simulation data of NUMT mtDNA coverage and perform pairwise comparisons between 320 windows representing 16,000 bp in mtDNA. The input data (list) contains the mtDNA coverage for each of 45 species.


#Length of NUMTs for each species. The information is stored in the file 'NUMT_length_species.txt'.

list=list(Shar,	Pdis,	Mmyo,	Pkuh,	Mmol,	Raeg,	Rfer,	Btau,	Oorc,	Sscr,	Cfer,	Csim,	Ecab,	Mjav,	Fcat,	Hhya,	Cfam,	Mput,	Zcal,	Ccri,	Eeur,	Sara,	Moch,	Mmus,	Rnor,	Hgla,	Svul,	Ocun,	Opri,	Tbel,	Mmur,	Ogar,	Sbol,	Ptro,	Hsap,	Mfas,	Dnov,	Tman,	Lafr,	Eedw,	Etel,	Mdom,	Vurs,	Glea,	Oana) 

library(reshape2)

# Generate the function to deal with circularity issue of mtDNA. Any simulated NUMT position which is larger than 16,000 bp will be deducted by 16,000

test<-function(val) {ifelse(val>16000, abs(16000 -val), val)}  

# Repeat the analysis 1,000 times

y<-1
for(y in 1) {
  
  print(y)

# Generate random mtDNA coordinates of NUMTs

for(a in 1:length(list)) {
  
  xc<-sample(1:16000,length(list[[a]]),replace=FALSE,prob=c(rep(1/16000,200),rep(1/16000,15800)))
  end<-xc+list[[a]]
  
  vec<-vector()

# Generate simulated mtDNA coordinates of NUMTs for each species (Loop through all 45 species) and obtain the new (simulated distribution)
  
  START<-1
  
  for (i in 1:length(list[[a]])) {
    val<-xc[i]:(end[i]-1)
    vec[START:(length(val)+START-1)]<-val
    START<-START+length(val)
  }
  
  test(end)
  
  cova<-test(vec)
  
  covb<-table(factor(cova,level=c(1:16000)))

# Arrange the coverage table into 320 windows (50 bp per window)

  d<-as.data.frame(covb)
  e<-matrix(d[,2],ncol=320)

# Calculate the median for each window
  
  for(z in 1:ncol(e))
    
  {
    median=median(e[,z])
    write.table(median,"median_each.txt", sep="\t",append=TRUE, row.names=FALSE, col.names=FALSE)
  }
  
}

# Arrange data format suitable for pairwise Mann-Whitney U test

  m<-read.table("median_each.txt")
  v<-matrix(m[,1],nrow=320)
  write.table(v,"v.txt")
  rownames(v)<-c(1:320)
  colnames(v)<-c(1:45)
  melt<-melt(v,id.vars='')
  windows<-factor(melt$Var2,levels=melt$Var1[1:320])

# Calculate pairwise Mann-Whitney U tests between each window

  re=pairwise.wilcox.test(melt$value,windows,p.adj="fdr",paired=TRUE)
  re1=c(re$p.value)
  re2=re1[!is.na(re1)]
  num<-length(re2[re2<0.01])
  num1<-length(re2[re2>=0.01])
  write.table(num,"simulation2.txt",sep="\t",append=TRUE,row.names = FALSE,col.names = FALSE)
  file.remove("median_each.txt")
 
}

