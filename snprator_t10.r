
##for backcompatibility here the old function working with tair10 data 

##  depending on where you put the data you might need to modfy the respective path
#load('~/data/additional_data/all_acc.rda')
#load('~/data/additional_data/anno_tair10_snps.rda')
#load('~/data/additional_data/tair10.rda')
#load('~/data/additional_data/snps_2029.rda')
## X.prefix is the prefix for the X.folder

snporator_t10<-function(gene,X.prefix='~/data/',data=1135,accessions=all,out.folder='~/data/snp_out') {
  if(!data%in%attributes(snporator_t10)$data) stop('Data need to be either 1135 or 2029')
  X.folder=paste(X.prefix,data,'_data/',sep='')
  
  setwd(X.folder)
  a1<-tair10[which(tair10[,5]==gene),2]
  
  a2<-tair10[which(tair10[,5]==gene),3]
  chr<-tair10[which(tair10[,5]==gene),1]
  
  
  D<-anno[which(anno[,1]==chr&anno[,2]>a1&anno[,2]<a2),]
  h<-1
  for ( r in 1:nrow(D)) {
    if(D[h,3]%in%SNPs$SNP==FALSE) { h=h+1} else {break}}
  
  load(paste('X_',data,'_',SNPs[which(SNPs$SNP==D[h,3]),5],'.rda',sep=''))
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
  

  
df<-D_[grep('NON_S',D_[,6]),]
  
get_isoforms <- function(snpeff, gene_id) {
  parts <- unlist(strsplit(snpeff, ","))
  isoforms <- unique(unlist(regmatches(parts, gregexpr(paste0(gene_id, "\\.[0-9]+"), parts))))
  isoforms
}


all_isoforms <- unique(unlist(lapply(df[[6]], get_isoforms, gene_id = gene)))

# Function: make one df per isoform with only NON_SYNONYMOUS for that isoform
make_isoform_df <- function(df, isoform) {
  keep_rows <- sapply(df[[6]], function(entry) grepl(isoform, entry))
  subdf <- df[keep_rows, ]
  
  subdf[[6]] <- sapply(subdf[[6]], function(entry) {
    parts <- unlist(strsplit(entry, ","))
    matches <- parts[grepl("NON_SYNONYMOUS", parts) & grepl(isoform, parts)]
    paste(matches, collapse = ",")
  })
  
  # Drop rows where no NON_SYNONYMOUS for that isoform was found
  subdf <- subdf[nchar(subdf[[6]]) > 0, ]
  
  subdf
}

# Generate a list of isoform-specific dfs
isoform_dfs <- lapply(all_isoforms, function(iso) make_isoform_df(df, iso))
names(isoform_dfs) <- all_isoforms

  ## now filter for the empty ones.....   
  
  A<-data.frame(SNP=c('length','total','non','syn','start','stop','low','moderate','high'),Number=c(a2-a1,nrow(D_),length(grep('NON_S',D_[,6])),(length(grep('SYN',D_[,6]))-length(grep('NON_S',D_[,6]))),length(grep('START',D_[,6])),length(grep('STOP',D_[,6])),length(grep('LOW',D_[,6])),length(grep('MODERATE',D_[,6])),length(grep('HIGH',D_[,6]))))
  if (dim(df)[1]>0) {
    V<-strsplit(df[,6],split="|",fixed=T)
    for ( i in 1: nrow(df)) {
      df[i,8]<-paste(unique(V[[i]][grep('NON_S',V[[i]])+3]),collapse=',')}
    if (dim(df)[1]>0) 
      df<-df[,-6]
    
    X2<-XX[,colnames(XX)%in%df[,3]]
    
    colnames(df)[7]<-'AA exchange'
    ## get accessions with alternative Allele 
    df$accessions<-NA
    for ( z in 1:nrow(df)) {
      df[z,8]<-paste(rownames(X2)[which(X2[,z]==1)],collapse=',')}
  }
  
  
  
setwd(out.folder)
  
  for (iso in names(isoform_dfs)) {
    outname <- paste0(iso, ".csv")
    write.csv(isoform_dfs[[iso]], outname, row.names = FALSE)
  }
  
  
snp_to_iso <- lapply(names(isoform_dfs), function(iso) {
  data.frame(SNP = isoform_dfs[[iso]]$SNP, isoform = iso, stringsAsFactors = FALSE)
})

snp_to_iso <- do.call(rbind, snp_to_iso)   

df_iso <- merge(df, aggregate(isoform ~ SNP, data = snp_to_iso, FUN = function(x) paste(unique(x), collapse = ",")),
                by = "SNP", all.x = TRUE)


  ####
  
  out1<-paste(gene,lol,'acc_snp_summary_t10.txt',sep='_')
  out2<-paste(gene,lol,'acc_non_syn_snps_t10.txt',sep='_')
  out3<-paste(gene,lol,'acc_snp_count_t10.txt',sep='_')

  write.table(A,file=out1,row.names=FALSE)
  write.table(df_iso,file=out2,row.names=FALSE)
  write.table(D2,file=out3,row.names=FALSE)
  
  cat('files',out1,out2,out3,'have been generated!','\n')
}

attr(snporator_t10,'data')<-c(1135,2029)

