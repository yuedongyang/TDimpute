rm(list = ls())
print('start...')

if(!library("survival", logical.return = TRUE)){
  install.packages("survival")
}
if(!library("plyr", logical.return = TRUE)){
  install.packages("plyr")
}
if(!library("cluster", logical.return = TRUE)){
  install.packages("cluster")
}
if(!library("aricode", logical.return = TRUE)){
  install.packages("aricode")
}

args<-commandArgs(T)
cancertype <- 'WT' 
method <- as.numeric(args[1])  
n_clusters<- 2 
missing_perc <- 0.5
sample_count <- 2
cat(cancertype, n_clusters, method,'\n')
datadir <- '/data' ## home directory

###########read clinical data
clin_file <- paste(datadir,"/data_clinical_patient.txt", sep="")
datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
clin_info <- datClin[c('PATIENT_ID','OS_DAYS','OS_STATUS')]
clin_info <- as.data.frame(clin_info)			
colnames(clin_info) <- c('bcr_patient_barcode','time','status')	
clin_info$status <- ifelse(clin_info$status == 'LIVING', 0,1)
clin <- clin_info[c('bcr_patient_barcode','time','status')]		
clin <- na.omit(clin)


###########read imputed data 
if(method==1){
  expr_file<-paste(datadir,"/filled_data/TDimpute_",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM_ref_1.csv",sep = "") 
}else if(method==2){
  expr_file<-paste(datadir,"/filled_data/TDimpute_self_",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM_ref_1.csv",sep = "") 
}else if(method==3){
  expr_file<-paste(datadir,"/filled_data/SVD_RNA_DNA_",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM_ref_1.csv",sep = "") 
}else if(method==4){
  expr_file<-paste(datadir,"/filled_data/TOBMI_",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM_ref_1.csv",sep = "") 
}else if(method==5){
  expr_file<-paste(datadir,"/filled_data/Lasso_RNA_DNA_",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM_ref_1.csv",sep = "") 
}else{
  expr_file<-paste("") 
}
print(expr_file)
datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
rnames<- datExpr[,1]
datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
datExpr <- na.omit(datExpr)
datExpr<-datExpr[,1:18155]  #only gene expression matrix
datExpr_original <- datExpr[!duplicated(row.names(datExpr)), ]  # remove duplicated samples
datExpr<-datExpr_original
datExpr[datExpr<0] <- 0

#merge expression matrix with clinical matrix
datExpr <- datExpr[substr(rownames(datExpr),18,19) == '01', ] # 11 and 06 samples are removed
rownames(datExpr) <- substr(rownames(datExpr), 1, 16)
datExpr <- as.data.frame(datExpr)
datExpr$bcr_patient_barcode <- row.names(datExpr)
merged_data <- merge(clin,datExpr,by.x="bcr_patient_barcode",by.y="bcr_patient_barcode")
merged_data <- merged_data[merged_data$time>0,]
exprSet<-merged_data
gene_names <- colnames(exprSet[4:ncol(exprSet)])
gene_names <- chartr("|",".",gene_names)
gene_names <- chartr("-",".",gene_names)
gene_names <- chartr("?","X",gene_names)
colnames(exprSet)<-c(colnames(exprSet[1:3]), gene_names)
cat('dim(exprSet):  ',dim(exprSet),'\n')	

# univariate Cox for cancer-related gene
mysurv <- Surv(as.numeric(exprSet$time), exprSet$status)
Unicox <- function(x){
  fml <- as.formula(paste0('mysurv~', x))
  gcox <- coxph(fml, exprSet)
  cox_sum <- summary(gcox)
  HR <- round(cox_sum$coefficients[,2],2)
  PValue <- round(cox_sum$coefficients[,5],4)
  CI <- paste0(round(cox_sum$conf.int[,3:4],2),collapse='-')
  Uni_cox <- data.frame('Characteristics' = x,
                        'Hazard.Ratio' = HR,
                        'CI95' = CI,
                        'P.value' = PValue)
  return(Uni_cox)
}
VarNames <- colnames(exprSet[,4:ncol(exprSet)]) #gene_names
Univar <- lapply(VarNames, Unicox)
Univar <- ldply(Univar, data.frame)
sort_genes <- Univar[order(Univar$P.value),]
sign_gene <- sort_genes[1:100,]	
sign_gene <- na.omit(sign_gene)
sign_gene_name<-as.character(sign_gene$Characteristics)
df_expr <- exprSet[,sign_gene_name]

##divided into 2 clusters by K-means
k<-n_clusters
set.seed(123)
km_res <- kmeans(df_expr, k, iter.max=100, nstart=200)
ss<-silhouette(km_res$cluster, dist(df_expr))
cat(mean(ss[, 3]), km_res$size,'\n')
yy <- Surv(as.numeric(exprSet$time), exprSet$status)
risk_level<-as.factor(ifelse(km_res$cluster>1,"High","Low"))
kms<-survfit(yy~risk_level,data=as.data.frame(exprSet))
kmdffexp=survdiff(yy~risk_level,data=as.data.frame(exprSet))
pValue=1-pchisq(kmdffexp$chisq, df=1)
cat('pValue result: ', pValue,'\n')

aa<-exprSet[,1:3]
yy <- Surv(as.numeric(aa$time), aa$status)
risk_level<-as.factor(ifelse(km_res$cluster>1,"Group 1","Group 2"))
kms<-survfit(yy~risk_level,data=as.data.frame(aa))
kmdffexp=survdiff(yy~risk_level,data=as.data.frame(aa))
pValue=format(1-pchisq(kmdffexp$chisq,df=1), digits = 3,  scientific = TRUE) 
palette=c("red","green")
ggsurvplot(kms, conf.int=F,xlab = "Time (days)",pval = paste("p = ",pValue,sep = ""),pval.coord = c(3000, 1),font.legend=14,
           font.x=14, font.y=14, palette=palette,  font.tickslab=14,legend.title ="", legend.labs=c("High risk", "Low risk"), title='TDimpute')

