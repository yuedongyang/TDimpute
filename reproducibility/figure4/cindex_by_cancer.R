###compare survival analyses for datasets imputed by different methods, Figure 4B in the paper

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
reduced_names = c('BRCA', "SARC", "LUSC", "BLCA", "KIRC", "LGG_", "PRAD", "LUAD", "LIHC", "SKCM",
                 "HNSC", "CESC", "COAD", "KIRP", "THCA", "STAD")
				 
args<-commandArgs(T)
p_cutoff<- 0.05
alpha_cox<-0 # ridge regression, alpha<-0.5 elastic net
method <- as.numeric(args[1])  # TDimpute:1, TDimpute_self:2, SVD:3, TOBMI:4, Lasso:5
datadir <- '/data' 


for(cancertype in reduced_names){
	##read clinical data 
	clin_file <- paste(datadir,"/clinical_cancers/test__nationwidechildrens.org_clinical_patient_",tolower(cancertype),".txt", sep="")
	datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
	clin_info <-datClin[c('bcr_patient_barcode','days_to_death','vital_status','days_to_last_followup')]
	clin_info <- as.data.frame(clin_info)
	clin_info$time <- ifelse(clin_info$vital_status == 'Alive', clin_info$days_to_last_followup,clin_info$days_to_death)
	clin_info$status <- ifelse(clin_info$vital_status == 'Alive', 0,1)
	clin <- clin_info[c('bcr_patient_barcode','time','status')]
			
    results<-matrix(nrow = 5, ncol=2) 
	pp_test<-array(0, dim=c(5, 5, 1000)) 
	
    perc<-1
    print(cancertype)
    for(missing_perc in c(0.1,0.3,0.5,0.7,0.9)){ 
        results_by_perc<-matrix(nrow = 5, ncol=2)         
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
			
			## 5-fold cross-validation
			fold_n<-5
			results_train<-matrix(nrow = 1, ncol=fold_n) 
			results_test<-matrix(nrow = 1, ncol=fold_n)
			set.seed(30)  
			folds<-createFolds(exprSet$bcr_patient_barcode,k=fold_n)				
			fold_count<-1
			risk_level_total<-c()
			test_data_total<-data.frame()
			for(i in 1:fold_n){				
				test_data<-exprSet[folds[[i]],]
				train_data<-exprSet[-folds[[i]],]								
				test_data_total<-rbind(test_data_total, test_data)
				cat('train size and test size:',nrow(train_data),nrow(test_data),'\n')

				rem <- function(x){
				r <-as.numeric(apply(x,2,function(i) sum(i>0)))
				  selected <- which(r > dim(x)[1]*0.5)
				  return(selected)
				}
				dataset_aa<-train_data[,4:ncol(train_data)]
				selected<- rem(dataset_aa)
				dataset_aa<- dataset_aa[,selected]
				train_data <- cbind(train_data[,1:3], dataset_aa)
				cat('dim(train_data): ',dim(train_data),'\n')
				
				# univariate Cox for cancer-related gene
				mysurv <- Surv(as.numeric(train_data$time), train_data$status)				
				Unicox <- function(x){
					flag<-0
					result<-tryCatch({
					  fml <- as.formula(paste0('mysurv~', x))
					  gcox <- coxph(fml, train_data)
					  cox_sum <- summary(gcox)
					  HR <- round(cox_sum$coefficients[,2],2)
					  PValue <- round(cox_sum$coefficients[,5],4)
					  CI <- paste0(round(cox_sum$conf.int[,3:4],2),collapse='-')
					  Uni_cox <- data.frame('Characteristics' = x,
											'Hazard.Ratio' = HR,
											'CI95' = CI,
											'P.value' = PValue)
					   flag<-1
					}, error   = function(e) { 
					  cat("error:",x,flag,"\n")
					}, finally = {
					   if(flag==0){
						  Uni_cox <- data.frame('Characteristics' = x,
												'Hazard.Ratio' = 0,
												'CI95' = 0,
												'P.value' = 1.0)			   
					   }
					   return(Uni_cox)
					})
					return(result)
				}					
				VarNames <- colnames(train_data[,4:ncol(train_data)])
				Univar <- lapply(VarNames, Unicox)
				Univar <- ldply(Univar, data.frame)
				sort_genes = Univar[order(Univar$P.value),]			
				sign_gene = sort_genes[sort_genes$P.value<p_cutoff, ]
				if(dim(sign_gene)[1]==0){
					top_n<-100
					sign_gene <- sort_genes[1:top_n, ]
				}
				cat('sign_gene size:',dim(sign_gene),'\n')
				sign_gene<-na.omit(sign_gene)
				sign_gene_name<-as.character(sign_gene$Characteristics)

				##filter by cancer-related gene
				x <- as.matrix(train_data[,sign_gene_name])
				y <- Surv(as.numeric(train_data$time), train_data$status)

				#registerDoParallel(cores=20)
				set.seed(1011) 
				cvfit<-cv.glmnet(x,y,family="cox", alpha=alpha_cox,nfold=10, maxit = 10000)  #alpha=1,lasso-cox
				#stopImplicitCluster()
				
				pp<-predict(cvfit, as.matrix(train_data[,sign_gene_name]), s=cvfit$lambda.min) #lambda.1se
				yy <- Surv(as.numeric(train_data$time), train_data$status)
				c_index<-survConcordance(formula = yy ~ pp)
				cat('c-index, train:',round(c_index$concordance,2),'\n')
				results_train[1,fold_count]<- round(c_index$concordance,3) 		
				
				pp<-predict(cvfit, as.matrix(test_data[,sign_gene_name]), s=cvfit$lambda.min) #lambda.1se
				yy <- Surv(as.numeric(test_data$time), test_data$status)
				c_index<-survConcordance(formula = yy ~ pp)
				cat('c-index, test:',round(c_index$concordance,2),'\n')			
				#results_test[1,fold_count]<- round(c_index$concordance,3) 
				if(is.na(c_index$concordance)){
					results_test[1,fold_count]<- 0.5
				}else{
					results_test[1,fold_count]<- round(c_index$concordance,3)
				}
				if(sort_genes[1,]$P.value<0){
					results_test[1,fold_count]<- 0.5
				}else{
					results_test[1,fold_count]<- round(c_index$concordance,3) 
				}
				fold_count <- fold_count + 1
			}
			cat('mean by fold:',round(rowMeans(results_train),3),'\n')
			cat('mean by fold:',round(mean(na.omit(results_test[1,])),3),'\n')
			results_by_perc[1,]<-c( round(rowMeans(results_train),3), round(mean(na.omit(results_test[1,])),3))
		}
	    cindex_file<-paste(datadir,"/cindex/method_",method,"_", cancertype, 5, missing_perc*100, ".RData",sep = "")
	    print(cindex_file)
	    save(results_by_perc, file=cindex_file)
	   
		results[perc,]<- colMeans(results_by_perc)
		print(results)
		perc<-perc+1
	}			
    print(results)
}

