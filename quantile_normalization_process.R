########## quantile normalization for DNA methylation from different source (e.g. TCGA and TARGET)
expr_file <- '/data0/zhoux/TARGET/DNA_WT.csv'
dat_target <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
aa <-data.matrix(dat_target[,2:dim(dat_target)[2]])
cnames<-colnames(dat_target)[2:dim(dat_target)[2]]
rnames<- dat_target[,1]
dat_target<-matrix(aa,nrow=dim(dat_target)[1],ncol=dim(dat_target)[2]-1,dimnames=list(rnames,cnames))


library("preprocessCore")
# save(target, file='/data0/zhoux/reference_distribution_DNA.RData')
load(file='/data0/zhoux/reference_distribution_DNA.RData')

WT.normalized <- normalize.quantiles.use.target(t(dat_target), target)
WT.normalized <- t(WT.normalized)
row.names(WT.normalized) <- row.names(dat_target)
colnames(WT.normalized) <- colnames(dat_target)
write.csv(WT.normalized, file='/data0/zhoux/TARGET/quantiles_DNA_WT_RSEM.csv', quote=F)
###############



######### quantile normalization for RNA-seq data from different source (e.g. TCGA and TARGET)
expr_file <- '/data0/zhoux/TARGET/UCSC_RNA_WT.csv'
dat_target <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
aa <-data.matrix(dat_target[,2:dim(dat_target)[2]])
cnames<-colnames(dat_target)[2:dim(dat_target)[2]]
rnames<- dat_target[,1]
dat_target<-matrix(aa,nrow=dim(dat_target)[1],ncol=dim(dat_target)[2]-1,dimnames=list(rnames,cnames))
var_target <- apply(dat_target, 2, var)>0
dat_target_filter <- dat_target[, var_target]

library("preprocessCore")
# save(target, file='reference_distribution_RNA.RData')
load(file='/data0/zhoux/reference_distribution_RNA.RData')

WT.normalized <- normalize.quantiles.use.target(t(dat_target_filter),target)
WT.normalized <- t(WT.normalized)
row.names(WT.normalized) <- row.names(dat_target_filter)
colnames(WT.normalized) <- colnames(dat_target_filter)
write.csv(WT.normalized, file='/data0/zhoux/TARGET/quantiles_RNA_WT_RSEM_ref_1.csv', quote=F)
####################