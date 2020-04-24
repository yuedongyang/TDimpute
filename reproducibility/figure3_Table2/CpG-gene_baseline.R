###compute baseline pearson correlation between RNAseq matrix and DNA methylation matrix
###run before CpG-gene_compare.R and CpG-gene_select.R
rm(list = ls())
print('start...')
library('psych')
options(stringsAsFactors=F)

args<-commandArgs(T)
cancertype<- args[1] #'BLCA' 
print(cancertype)
datadir <- '/data'
expr_file<-paste(datadir,"/bootstrap_cancer_V1/",cancertype,"10_1.csv",sep = "") 
print(expr_file)
shuffle_cancer <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
aa <-data.matrix(shuffle_cancer[,2:dim(shuffle_cancer)[2]])
cnames<-colnames(shuffle_cancer)[2:dim(shuffle_cancer)[2]]
rnames<- shuffle_cancer[,1]
shuffle_cancer<-matrix(aa,nrow=dim(shuffle_cancer)[1],ncol=dim(shuffle_cancer)[2]-1,dimnames=list(rnames,cnames))
rnames<- row.names(shuffle_cancer)
shuffle_cancer <- shuffle_cancer[!duplicated(rnames), ] 

#RNA_size: 19027
#DNA_size: 27717
baseline_coeff<-cor( x=shuffle_cancer[,1:19027], y=shuffle_cancer[,19028:ncol(shuffle_cancer)] ) #compute pearson correlation
baseline_coeff <- na.omit(baseline_coeff) ###drop the NA gene-CpG pair
cat('dim(baseline_coeff): ', dim(baseline_coeff), '\n')
max_value<-vector(mode = "numeric", length = nrow(baseline_coeff))
max_CpG<-vector(mode = "character", length = nrow(baseline_coeff))
col_names <- colnames(baseline_coeff)
for(i in c(1:nrow(baseline_coeff))){
	max_value[i] <- max( baseline_coeff[i,]**2 )  #pick the strongest pair
	max_CpG[i] <- col_names[which.max(baseline_coeff[i,]**2)]
}
max_pair<-data.frame('gene'=row.names(baseline_coeff), 'max_pcc'=max_value, 'CpG'=max_CpG) 

max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
for(i in c(1:nrow(max_pair))){
	res<-cor.test( shuffle_cancer[, max_pair[i,]$gene], shuffle_cancer[, max_pair[i,]$CpG] )
	max_pvalue[i] <- res$p.value
}
max_pvalue <- p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue)) #FDR-q
max_pair$adj_pvalue<-max_pvalue

coeff_file<-paste(datadir,"/coeff_matrix/coeff_matrix_base",cancertype,"_gene_CPG_pair.RData",sep = "")
print(coeff_file)
save(max_pair, file=coeff_file)


