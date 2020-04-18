### extracting SNP information for 1001 genomes data genewise, tair10 annotation  
## all data needed are available at go.uniwue.de/snpextractor 
## can be run with resequencing data (data=1135, default) or with imputed 2029 data (data=2029)

## will change script to use non_imputed data as well in the future

## the gene needs to be specified with standard AGI code (gene<-'AT4G04740')


##  depending on where you put the data you might need to modfy the respective path to load the required data 
#load('data/additional_data/all_acc.rda')
#load('data/additional_data/anno_ara11_snps.rda')
#load('data/additional_data/ara11.rda')
#load('data/additional_data/snps_2029.rda')
## X.prefix is the prefix for the X.folder


wd<-getwd()

snporator<-function(gene,X.prefix='~/data/',data=1135,accessions=all,out.folder=wd) {
if(!data%in%attributes(snporator)$data) stop('Data need to be either 1135 or 2029')
X.folder=paste(X.prefix,data,'_data/',sep='')

setwd(X.folder)



D<-Z_[which(Z_$gene==gene),]
if(any(D[,3]%in%SNPs$SNP)==T) {
h<-1
 for ( r in 1:nrow(D)) {
if(D[h,3]%in%SNPs$SNP==FALSE) { h=h+1} else {break}}

load(paste('X_',data,'_',SNPs[which(SNPs$SNP==D[h,3]),5],'.rda',sep=''))
XX<-as.matrix(X[as.numeric(rownames(X))%in%accessions,colnames(X)%in%D$SNP])
colnames(XX)=D$SNP
rm(X)
lol<-nrow(XX)
j<-nrow(D)
 for ( r in 1:nrow(D)) {
if(D[j,3]%in%SNPs$SNP==FALSE) { j=j-1} else {break}}


 if(SNPs[which(SNPs$SNP==D[j,3]),5]!=SNPs[which(SNPs$SNP==D[h,3]),5]) {
 
load(paste('X_',data,'_',SNPs[which(SNPs$SNP==D[nrow(D),3]),5],'.rda',sep=''))
XX<-cbind(XX,X[as.numeric(rownames(X))%in%accessions,colnames(X)%in%D$SNP])
rm(X) }

D_<-subset(D,D$SNP%in%colnames(XX))

D_$count<-apply(XX,2,sum)

 A1<-table(D_$variant)
 ## output start,stop and missence, modify if you need something else 
 Dns<-D_[grep('missense_variant',D_$variant),]
 D_2<-D_[grep('start',D_$variant),]
 if(nrow(D_2)>0) {
         Dns<-rbind(Dns,D_2)}
 
 D_3<-D_[grep('stop',D_$variant),]
 if(nrow(D_3)>0) {
         Dns<-rbind(Dns,D_3)}
 
X2<-as.matrix(XX[,colnames(XX)%in%Dns[,3]])
## get accessions with alternative Allele 
if(nrow(Dns)>0) {
Dns$accessions<-NA
 for ( z in 1:nrow(Dns)) {
 Dns[z,11]<-paste(rownames(X2)[which(X2[,z]==1)],collapse=',')}}



out1<-paste(gene,lol,'snp_summary_a11.txt',sep='_')
out2<-paste(gene,lol,'impact_snps_peracc_a11.txt',sep='_')
out3<-paste(gene,lol,'snp_count_a11.txt',sep='_')
setwd(out.folder)
write.table(A1,file=out1,row.names=FALSE)
write.table(Dns,file=out2,row.names=FALSE)
write.table(D_,file=out3,row.names=FALSE)

cat('files',out1,out2,out3,'have been generated!','\n')
} else {cat('no snps for the gene',gene,'have been detected')}

}

attr(snporator,'data')<-c(1135,2029)
