###compute baseline for prognosis-related genes

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
if(!library("PRROC", logical.return = TRUE)){
	install.packages("PRROC")
}

reduced_names = c('BRCA', "SARC", "LUSC", "BLCA", "KIRC", "LGG_", "PRAD", "LUAD", "LIHC", "SKCM",
                 "HNSC", "CESC", "COAD", "KIRP", "THCA", "STAD")
datadir <- '/data' 

for(cancertype in reduced_names){
    print(cancertype)
	##read imputed data 
	clin_file <- paste(datadir,"/clinical_cancers/test__nationwidechildrens.org_clinical_patient_",tolower(cancertype),".txt", sep="")
	datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
	clin_info <-datClin[c('bcr_patient_barcode','days_to_death','vital_status','days_to_last_followup')]
	clin_info <- as.data.frame(clin_info)
	clin_info$time <- ifelse(clin_info$vital_status == 'Alive', clin_info$days_to_last_followup,clin_info$days_to_death)
	clin_info$status <- ifelse(clin_info$vital_status == 'Alive', 0,1)
	clin <- clin_info[c('bcr_patient_barcode','time','status')]

    ######################################
    expr_file<-paste(datadir, "/bootstrap_cancer_V1/",cancertype,"10_1.csv",sep = "")
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

    rem <- function(x){
      r <-as.numeric(apply(x,2,function(i) sum(i>0)))
      selected <- which(r > dim(x)[1]*0.5)
      return(selected)
    }
    selected<- rem(datExpr)
    datExpr<- datExpr[,selected]

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

    ################coxph model########
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
    VarNames <- gene_names
    Univar <- lapply(VarNames, Unicox)
    Univar <- ldply(Univar, data.frame)
    Univar[,5]=p.adjust(Univar$P.value, method ="fdr", n=dim(Univar)[1])
    colnames(Univar)<-c("Characteristics", "Hazard.Ratio", "CI95", "P.value", "adj.p")

	sort_genes = Univar[order(Univar$P.value),]
    baseline_gene_file<-paste(datadir,'/cox_results/',cancertype,'_baseline_gene_list.csv',sep = "")
    cat('dim(sort_genes): ',dim(sort_genes),'\n')
    write.csv(sort_genes, file=baseline_gene_file , quote=F)
}

