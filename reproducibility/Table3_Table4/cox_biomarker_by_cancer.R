### compare the prognosis-related gene lists from different imputed datasets, Table 4 in the paper

rm(list = ls())
print('start...')

if(!library("survival", logical.return = TRUE)){
	install.packages("survival")
}
if(!library("plyr", logical.return = TRUE)){
	install.packages("plyr")
}
if(!library("caret", logical.return = TRUE)){
	install.packages("caret")
}
if(!library("glmnet", logical.return = TRUE)){
	install.packages("glmnet")
}

args<-commandArgs(T)
cancertype<- args[1] #e.g. 'BLCA'
prognostic <- args[2] #e.g. urothelial_BLCA
method <- as.numeric(args[3]) # TDimpute:1, TDimpute_self:2, SVD:3, TOBMI:4, Lasso:5
datadir <- '/data' 

###########read clinical data
clin_file <- paste(datadir,"/clinical_cancers/test__nationwidechildrens.org_clinical_patient_",tolower(cancertype),".txt", sep="")
datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
clin_info <-datClin[c('bcr_patient_barcode','days_to_death','vital_status','days_to_last_followup')]
clin_info <- as.data.frame(clin_info)
clin_info$time <- ifelse(clin_info$vital_status == 'Alive', clin_info$days_to_last_followup,clin_info$days_to_death)
clin_info$status <- ifelse(clin_info$vital_status == 'Alive', 0,1)
clin <- clin_info[c('bcr_patient_barcode','time','status')]


###########read original data 
expr_file<-paste(datadir,"/bootstrap_cancer_V1/",cancertype,"10_1.csv",sep = "") 
datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
rnames<- datExpr[,1]
datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
datExpr<-datExpr[,1:19027]#only gene expression matrix
datExpr_original <- datExpr[!duplicated(row.names(datExpr)), ]  # remove duplicated samples
datExpr<-datExpr_original
datExpr[datExpr<0] <- 0

datExpr_tumor <- datExpr[substr(rownames(datExpr),14,15) == '01'| substr(rownames(datExpr),14,15) == '06', ] # 11 samples are removed
rownames(datExpr_tumor) <- substr(rownames(datExpr_tumor), 1, 15)  #replicated sample like: 'TCGA-BK-A26L-01_Rep245'
datExpr_tumor<-datExpr_tumor[!duplicated(rownames(datExpr_tumor)),]
rownames(datExpr_tumor) <- substr(rownames(datExpr_tumor), 1, 12)
datExpr_tumor <- as.data.frame(datExpr_tumor)
row.names(datExpr_tumor) <- chartr(".","-",row.names(datExpr_tumor))
datExpr_tumor$bcr_patient_barcode <- row.names(datExpr_tumor)
merged_data<-merge(clin,datExpr_tumor,by.x="bcr_patient_barcode",by.y="bcr_patient_barcode")
merged_data <- merged_data[merged_data$time>0,]
exprSet<-merged_data
gene_names <- colnames(exprSet[4:ncol(exprSet)])
gene_names <- chartr("|",".",gene_names)
gene_names <- chartr("-",".",gene_names)
gene_names <- chartr("?","X",gene_names)
colnames(exprSet)<-c(colnames(exprSet[1:3]), gene_names)
cat('dim(exprSet):  ',dim(exprSet),'\n')

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

baseline_gene_file<-paste(datadir,'/cox_results/',cancertype,'_baseline_gene_including_06.csv',sep = "")
baseline_gene <- read.table(baseline_gene_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)

VarNames <- baseline_gene$Characteristics  ###gene_names
Univar <- lapply(VarNames, Unicox)
Univar <- ldply(Univar, data.frame)
Univar[,5]=p.adjust(Univar$P.value, method ="fdr", n=dim(Univar)[1])
colnames(Univar)<-c("Characteristics", "Hazard.Ratio", "CI95", "P.value", "adj.p")
sort_genes_base = Univar[order(Univar$P.value),]

gene_baseline_file<-paste(datadir,'/baseline_prognostic_biomarker/prognostic_',prognostic,'.csv',sep = "")
gene_baseline <- read.table(gene_baseline_file,sep = ",",header=1,row.names=1, stringsAsFactors=FALSE,check.names=F)
gene_baseline <- gene_baseline[,1]

top_n<-100
test_gene = sort_genes_base[1:top_n, ]
top100_baselinegene = c()
for(i in 1:top_n){
	top100_baselinegene[i]<-unlist(strsplit(as.character(test_gene[i,'Characteristics']),"\\."))[1]
}
True_comm_gene_1<-intersect(gene_baseline, top100_baselinegene)
cat('True_comm_gene_1: ',length(True_comm_gene_1), '\n')
	
top_n<-200
test_gene = sort_genes_base[1:top_n, ]
top200_baselinegene = c()
for(i in 1:top_n){
	top200_baselinegene[i]<-unlist(strsplit(as.character(test_gene[i,'Characteristics']),"\\."))[1]
}
True_comm_gene_2<-intersect(gene_baseline, top200_baselinegene)
cat('True_comm_gene_2: ',length(True_comm_gene_2), '\n')	
######################################################################################


clin_file <- paste(datadir,"/clinical_cancers/test__nationwidechildrens.org_clinical_patient_",tolower(cancertype),".txt", sep="")
datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
clin_info <-datClin[c('bcr_patient_barcode','days_to_death','vital_status','days_to_last_followup')]
clin_info <- as.data.frame(clin_info)
clin_info$time <- ifelse(clin_info$vital_status == 'Alive', clin_info$days_to_last_followup,clin_info$days_to_death)
clin_info$status <- ifelse(clin_info$vital_status == 'Alive', 0,1)
clin <- clin_info[c('bcr_patient_barcode','time','status')]

results=matrix(nrow = 5, ncol=6) # ncol=4 indicate top_100 or 200, comm_gene with True_comm_gene 
perc<-1
print(cancertype)
for(missing_perc in c(0.1,0.3,0.5,0.7,0.9)){  
	results_by_perc=matrix(nrow = 5, ncol=6)         
	for(sample_count in c(1,2,3,4,5)){
		##read imputed data 
		if(method==1){
		  expr_file<-paste(datadir,"/filled_data/TDimpute_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
		}else if(method==2){
		  expr_file<-paste(datadir,"/filled_data/TDimpute_self_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
		}else if(method==3){
		  expr_file<-paste(datadir,"/filled_data/SVD_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
		}else if(method==4){
		  expr_file<-paste(datadir,"/filled_data/TOBMI_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
		}else if(method==5){
		  expr_file<-paste(datadir,"/filled_data/Lasso_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
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
		datExpr<-datExpr[,1:19027] #only gene expression matrix
		datExpr_original <- datExpr[!duplicated(row.names(datExpr)), ]  # remove duplicated samples
		datExpr<-datExpr_original
		datExpr[datExpr<0] <- 0
		datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:ceiling(dim(datExpr)[2]*0.9)]]

		datExpr_tumor <- datExpr[substr(rownames(datExpr),14,15) == '01' | substr(rownames(datExpr),14,15) == '06', ] # 11 samples are removed
		rownames(datExpr_tumor) <- substr(rownames(datExpr_tumor), 1, 15)  #replicated sample like: 'TCGA-BK-A26L-01_Rep245'
		datExpr_tumor<-datExpr_tumor[!duplicated(rownames(datExpr_tumor)),]
		rownames(datExpr_tumor) <- substr(rownames(datExpr_tumor), 1, 12)
		datExpr_tumor <- as.data.frame(datExpr_tumor)
		row.names(datExpr_tumor) <- chartr(".","-",row.names(datExpr_tumor))
		datExpr_tumor$bcr_patient_barcode <- row.names(datExpr_tumor)
		merged_data<-merge(clin,datExpr_tumor,by.x="bcr_patient_barcode",by.y="bcr_patient_barcode")
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

		VarNames <- baseline_gene$Characteristics  ###gene_names
		VarNames <- intersect(colnames(exprSet)[4:ncol(exprSet)], VarNames)
		Univar <- lapply(VarNames, Unicox)
		Univar <- ldply(Univar, data.frame)
		Univar[,5] <- p.adjust(Univar$P.value, method ="fdr", n=dim(Univar)[1])
		colnames(Univar)<-c("Characteristics", "Hazard.Ratio", "CI95", "P.value", "adj.p")		
		sort_genes = Univar[order(Univar$P.value),]
		
		top_n<-100
		test_gene_1 = sort_genes[1:top_n, ]			
		genelist1 = c()
		for(i in 1:top_n){
			genelist1[i]<-unlist(strsplit(as.character(test_gene_1[i,'Characteristics']),"\\."))[1]
		}
		comm_gene1<-intersect(gene_baseline, genelist1)
		comm_gene1_true<-intersect(True_comm_gene_1, genelist1)
		cat('comm_gene1-> ','sample_count:',sample_count,length(comm_gene1),length(comm_gene1_true), length( intersect(top100_baselinegene, genelist1) ), '\n')			
		
		top_n<-200
		test_gene_1 = sort_genes[1:top_n, ]			
		genelist2 = c()
		for(i in 1:top_n){
			genelist2[i]<-unlist(strsplit(as.character(test_gene_1[i,'Characteristics']),"\\."))[1]
		}
		comm_gene2<-intersect(gene_baseline, genelist2)
		comm_gene2_true<-intersect(True_comm_gene_2, genelist2)
		cat('comm_gene2-> ','sample_count:',sample_count,length(comm_gene2), length(comm_gene2_true), length( intersect(top200_baselinegene, genelist2) ),'\n')		
		
		results_by_perc[sample_count,]<-c( length(comm_gene1), length(comm_gene2), length(comm_gene1_true), length(comm_gene2_true), length(intersect(top100_baselinegene, genelist1)), length(intersect(top200_baselinegene, genelist2)) )
	}
	coeff_file<-paste(datadir,"/cox_filter/biomarker_method_", method,"_", 5, cancertype, missing_perc*100, ".RData",sep = "")
	print(coeff_file)
	save(results_by_perc, file=coeff_file)
	
   results[perc,]<- round(apply(results_by_perc,2,mean),2)
   print(results)
   perc<-perc+1
}
print(results)



