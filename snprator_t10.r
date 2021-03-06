
##for backcompatibility here the old function working with tair10 data 

##  depending on where you put the data you might need to modfy the respective path
load('~/data/additional_data/all_acc.rda')
load('~/data/additional_data/anno_tair10_snps.rda')
load('~/data/additional_data/tair10.rda')
load('~/data/additional_data/snps_2029.rda')
## X.prefix is the prefix for the X.folder

snporator_t10<-function(gene,X.prefix='~/data/',data=1135,accessions=all,out.folder=wd) {
  if(!data%in%attributes(snporator)$data) stop('Data need to be either 1135 or 2029')
  X.folder=paste(X.prefix,data,'_data/',sep='')
  
  setwd(X.folder)
  a1<-tair10[which(tair10[,5]==gene),2]
  
  a2<-tair10[which(tair10[,5]==gene),3]
  chr<-tair10[which(tair10[,5]==gene),1]
  
  
  D<-anno[which(anno[,1]==chr&anno[,2]>a1&anno[,2]<a2),]
  h<-1
  for ( r in 1:nrow(D)) {
    if(D[h,3]%in%SNPs$SNP==FALSE) { h=h+1} else {break}}
  
  load(paste('X_',data,'_',SNPs[which(SNPs$SNP==D[h,3]),5],'.rda.gz',sep=''))
  XX<-X[as.numeric(rownames(X))%in%accessions,colnames(X)%in%D$SNP]
  rm(X)
  lol<-nrow(XX)
  j<-nrow(D)
  for ( r in 1:nrow(D)) {
    if(D[j,3]%in%SNPs$SNP==FALSE) { j=j-1} else {break}}
  
  
  if(SNPs[which(SNPs$SNP==D[j,3]),5]!=SNPs[which(SNPs$SNP==D[h,3]),5]) {
    
    load(paste('X_',data,'_',SNPs[which(SNPs$SNP==D[nrow(D),3]),5],'.rda.gz',sep=''))
    XX<-cbind(XX,X[as.numeric(rownames(X))%in%accessions,colnames(X)%in%D$SNP])
    rm(X) }
  
  D_<-subset(D,D$SNP%in%colnames(XX))
  
  D_$count<-apply(XX,2,sum)
  
  D2<-D_[,-6]
  
  Dns<-D_[grep('NON_S',D_[,6]),]
  
  A<-data.frame(SNP=c('length','total','non','syn','start','stop','low','moderate','high'),Number=c(a2-a1,nrow(D_),length(grep('NON_S',D_[,6])),(length(grep('SYN',D_[,6]))-length(grep('NON_S',D_[,6]))),length(grep('START',D_[,6])),length(grep('STOP',D_[,6])),length(grep('LOW',D_[,6])),length(grep('MODERATE',D_[,6])),length(grep('HIGH',D_[,6]))))
  if (dim(Dns)[1]>0) {
    V<-strsplit(Dns[,6],split="|",fixed=T)
    for ( i in 1: nrow(Dns)) {
      Dns[i,8]<-paste(unique(V[[i]][grep('NON_S',V[[i]])+3]),collapse=',')}
    if (dim(Dns)[1]>0) 
      Dns<-Dns[,-6]
    
    X2<-XX[,colnames(XX)%in%Dns[,3]]
    
    colnames(Dns)[7]<-'AA exchange'
    ## get accessions with alternative Allele 
    Dns$accessions<-NA
    for ( z in 1:nrow(Dns)) {
      Dns[z,8]<-paste(rownames(X2)[which(X2[,z]==1)],collapse=',')}
  }
  
  
  out1<-paste(gene,lol,'acc_snp_summary_t10.txt',sep='_')
  out2<-paste(gene,lol,'acc_non_syn_snps_t10.txt',sep='_')
  out3<-paste(gene,lol,'acc_snp_count_t10.txt',sep='_')
  setwd(out.folder)
  write.table(A,file=out1,row.names=FALSE)
  write.table(Dns,file=out2,row.names=FALSE)
  write.table(D2,file=out3,row.names=FALSE)
  
  cat('files',out1,out2,out3,'have been generated!','\n')
}

attr(snporator_t10,'data')<-c(1135,2029)

