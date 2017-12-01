## Permutation test to assess for enrichment of GPCRs in genes with monoallelic expression

# setwd("~/Desktop/Permutation_testing/Allelic_expression/") # Please set appropriate working directory

aexp = read.table('../data/allelicexp_genes.txt',header=T)
background<-as.vector(unlist(aexp$Ensemble_ID))

#Genes with mono-allelic expression
mae<-as.vector(unlist(aexp[which(aexp$MAE.1_BAE.0 == 1),1]))

#GPCR drug targets with mono-allelic expression
recep<-as.vector(unlist(read.table('GPCR_drugtargets.txt')))
recep.aexp = intersect(recep,background) # GPCR drug targets with allele-specific expression data
recep.mae<-intersect(recep.aexp,mae)

# Randomizations
random<-matrix(ncol=100000,nrow=length(recep.aexp))
for (i in 1:100000){
  sam<-sample(background, size=length(recep.aexp),replace=FALSE,prob=NULL)
  random[,i]<-sam
}

result<-vector(length=100000)
for(i in 1:ncol(random)){
  result[i]<-length(intersect(random[,i],mae))
}

l<-length(recep.mae)
m<-mean(result)
s<-sd(result)
z<-(length(recep.mae)-mean(result))/sd(result)
p<-length(which(result>=length(recep.mae)))/length(result)

res<-matrix(ncol=2, nrow=5)
res[,1]<-c("GPCRs with monoallelic expression","Mean of random expectation", "SD of random expectation", "Z-score","P-value")
res[,2]<-c(l,m,s,z,p)
write.table(res,'GPCRaexp_permut.txt',sep="\t",quote=F,row.names=F,col.names=F)
# if P-value is 0, then it should be reported as P<1E-5
