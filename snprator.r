### extracting SNP information for 1001 genomes data genewise, tair10 annotation 




#could change script to use non_imputed sata as well .
#gene<-'AT4G04740'
snporator<-function(gene,X.folder='/storage/evolgen/data/1135_data') {

out.folder<-getwd()

load('/storage/evolgen/data/anno_1135.rda')
load('/storage/evolgen/data/tair10.rda')
load('/storage/evolgen/data/snps_2029.rda')

setwd(X.folder)
a1<-tair10[which(tair10[,5]==gene),2]

a2<-tair10[which(tair10[,5]==gene),3]
chr<-tair10[which(tair10[,5]==gene),1]


D<-anno[which(anno[,1]==chr&anno[,2]>a1&anno[,2]<a2),]
h<-1
 for ( r in 1:nrow(D)) {
if(D[h,3]%in%SNPs$SNP==FALSE) { h=h+1} else {break}}

load(paste('X_1135_',SNPs[which(SNPs$SNP==D[h,3]),5],'.rda.gz',sep=''))
XX<-X[,colnames(X)%in%D$SNP]
rm(X)
j<-nrow(D)
 for ( r in 1:nrow(D)) {
if(D[j,3]%in%SNPs$SNP==FALSE) { j=j-1} else {break}}


 if(SNPs[which(SNPs$SNP==D[j,3]),5]!=SNPs[which(SNPs$SNP==D[h,3]),5]) {
 
load(paste('X_1135_',SNPs[which(SNPs$SNP==D[nrow(D),3]),5],'.rda.gz',sep=''))
XX<-cbind(XX,X[,colnames(X)%in%D$SNP])
rm(X) }

D_<-subset(D,D$SNP%in%colnames(XX))

D_$count<-apply(XX,2,sum)

D2<-D_[,-6]

 Dns<-D_[grep('NON_S',D_[,6]),]
 
A<-data.frame(SNP=c('length','total','non','syn','start','stop','low','moderate','high'),Number=c(a2-a1,nrow(D_),length(grep('NON_S',D_[,6])),(length(grep('SYN',D_[,6]))-length(grep('NON_S',D_[,6]))),length(grep('START',D_[,6])),length(grep('STOP',D_[,6])),length(grep('LOW',D_[,6])),length(grep('MODERATE',D_[,6])),length(grep('HIGH',D_[,6]))))
 
V<-strsplit(Dns[,6],split="|",fixed=T)
for ( i in 1: nrow(Dns)) {
 Dns[i,8]<-paste(unique(V[[i]][grep('NON_S',V[[i]])+3]),collapse=',')}

Dns<-Dns[,-6]

X2<-XX[,colnames(XX)%in%Dns[,3]]

colnames(Dns)[7]<-'AA exchange'
## get accessions with alternative Allele 
Dns$accessions<-NA
 for ( z in 1:nrow(Dns)) {
 Dns[z,8]<-paste(rownames(X2)[which(X2[,z]==1)],collapse=',')}



out1<-paste(gene,'snp_summary.txt',sep='_')
out2<-paste(gene,'non_syn_snps.txt',sep='_')
out3<-paste(gene,'snp_count.txt',sep='_')
setwd(out.folder)
write.table(A,file=out1,row.names=FALSE)
write.table(Dns,file=out2,row.names=FALSE)
write.table(D2,file=out3,row.names=FALSE)

cat('files',out1,out2,out3,'have been generated!','\n')
}

