## Permutation test to assess significance of overlap between known functional sites & putative functional sites

# setwd("~/Desktop/Permutation_testing/MV_permutation/") # Please set appropriate working directory
mut = read.table('../data/GPCR_Missense.txt',header=T) # all mutations
putat = as.numeric(unlist(mut[which(mut$func_putative == 'TRUE'),1]))
known = as.numeric(unlist(mut[which(mut$func_known == 'TRUE'),1]))
known.putat = intersect(known,putat)


background<-as.numeric(unlist(comp$Sno))

# Randomizations
random<-matrix(ncol=100000,nrow=length(known))
for (i in 1:100000){
  sam<-sample(background, size=length(known),replace=TRUE,prob=NULL)
  random[,i]<-sam
}

result<-vector(length=100000)
for(i in 1:ncol(random)){
  result[i]<-length(intersect(random[,i],putat))
}

l<-length(known.putat)
m<-mean(result)
s<-sd(result)
z<-(length(known.putat)-mean(result))/sd(result)
p<-length(which(result>=length(known.putat)))/length(result)


res<-matrix(ncol=2, nrow=5)
res[,1]<-c("Overlap between known and putative functional sites",
           "Mean of random expectation", "SD of random expectation", "Z-score","P-value")
res[,2]<-c(l,m,s,z,p)
write.table(res,'knownputat_permut.txt',sep="\t",quote=F,row.names=F,col.names=F)
