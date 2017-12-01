## Permutation test to assess significance of enrichment of known functional sites among sites with missense variants

# setwd("~/Desktop/Permutation_testing/MV_permutation/") # Please set appropriate working directory
sites = read.table('All_GPCRsites.txt',header=T) # ALL GPCR sites
table(sites$func_known)

mut = read.table('../data/GPCR_Missense.txt',header=T)
mut.func = as.numeric(unlist(mut[which(mut$func_known == "TRUE"),1]))
mut.ids = as.numeric(unlist(mut$Sno))

sites.func = as.numeric(unlist(sites[which(sites$func_known == 'TRUE'),]))

background<-as.numeric(unlist(sites$Sno))

# Randomizations

random<-matrix(ncol=100000,nrow=length(mut.ids))
for (i in 1:100000){
  sam<-sample(background, size=length(mut.ids),replace=TRUE,prob=NULL)
  random[,i]<-sam
}

result<-vector(length=100000)
for(i in 1:ncol(random)){
  result[i]<-length(intersect(random[,i],sites.func))
  print(i)
}

l<-length(mut.func)
m<-mean(result)
s<-sd(result)
z<-(length(mut.func)-mean(result))/sd(result)
p<-length(which(result>=length(mut.func)))/length(result)
#p<-length(which(result<=length(hr.ess)))/length(result)

res<-matrix(ncol=2, nrow=5)
res[,1]<-c("MVs at known functional sites","Mean of random expectation", "SD of random expectation", "Z-score","P-value")
res[,2]<-c(l,m,s,z,p)
write.table(res,'MVfuncknown_permut.txt',sep="\t",quote=F,row.names=F,col.names=F)
