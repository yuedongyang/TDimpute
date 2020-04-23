### compare the prognosis-related gene lists from different imputed datasets, Table 3 in the paper

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

args<-commandArgs(T)
p_cutoff<- 0.05
method <- as.numeric(args[1]) # TDimpute:1, TDimpute_self:2, SVD:3, TOBMI:4, Lasso:5
cat(p_cutoff, method,'\n')

for(cancertype in reduced_names){
	##read clinical data 
	clin_file <- paste(datadir,"/clinical_cancers/test__nationwidechildrens.org_clinical_patient_",tolower(cancertype),".txt", sep="")
	datClin <- read.delim(clin_file,sep = "\t",stringsAsFactors=FALSE,check.names=F)
	clin_info <-datClin[c('bcr_patient_barcode','days_to_death','vital_status','days_to_last_followup')]
	clin_info <- as.data.frame(clin_info)
	clin_info$time <- ifelse(clin_info$vital_status == 'Alive', clin_info$days_to_last_followup,clin_info$days_to_death)
	clin_info$status <- ifelse(clin_info$vital_status == 'Alive', 0,1)
	clin <- clin_info[c('bcr_patient_barcode','time','status')]
	
    results=matrix(nrow = 5, ncol=5) 
    perc<-1
    print(cancertype)
    for(missing_perc in c(0.1,0.3,0.5,0.7,0.9)){
        results_by_perc=matrix(nrow = 5, ncol=5)         
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

			baseline_gene_file<-paste(datadir,'/cox_results/',cancertype,'_baseline_gene_list.csv',sep = "")
            baseline_gene <- read.table(baseline_gene_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
            VarNames <- baseline_gene$Characteristics  ###gene_names
            Univar <- lapply(VarNames, Unicox)
            Univar <- ldply(Univar, data.frame)
            Univar[,5] <- p.adjust(Univar$P.value, method ="fdr", n=dim(Univar)[1])
            colnames(Univar)<-c("Characteristics", "Hazard.Ratio", "CI95", "P.value", "adj.p")
           
            sort_genes <- Univar[order(Univar$P.value),]
            base_cut <- baseline_gene
            pred_cut <- sort_genes
            true_pos <- base_cut$Characteristics[base_cut$P.value<p_cutoff]   
            true_neg <- base_cut$Characteristics[base_cut$P.value>=p_cutoff]  
            pred_pos <- pred_cut$Characteristics[pred_cut$P.value<p_cutoff]  
            pred_neg <- pred_cut$Characteristics[pred_cut$P.value>=p_cutoff]   
            pred_pos <- sub("_y","",pred_pos)
            pred_neg <- sub("_y","",pred_neg)  
                                   
            TP <- length(intersect(true_pos,pred_pos))
            FP <- length(pred_pos)-TP
            FN <- length(intersect(pred_neg,true_pos))

            PPV <- TP/(TP+FP)
            TPR <- TP/(TP+FN)
            Fscore <- (2*PPV*TPR)/(PPV+TPR)
            cat('true_pos, true_neg, pred_pos, pred_neg:',length(true_pos), length(true_neg), length(pred_pos), length(pred_neg), '\n')
            cat('TP, FP, FN:',TP,FP,FN, '\n')        
            cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')                                   
                       
            #### ROC curve- PR curve                       
            base_cut <- baseline_gene
			base_cut$P.value <- 1 - base_cut$P.value
			p_cutoff_new <- 1- p_cutoff
            base_cut$label[base_cut$P.value>p_cutoff_new] <- 1
            base_cut$label[base_cut$P.value<=p_cutoff_new] <- 0
            rownames(pred_cut)<-pred_cut$Characteristics
            rownames(base_cut)<-base_cut$Characteristics
            base_cut$pred_score <- 1-pred_cut[base_cut$Characteristics,'P.value']

            fg <- base_cut$pred_score[base_cut$label == 1]
            bg <- base_cut$pred_score[base_cut$label == 0]
            cat('..........before NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
			fg<-na.omit(fg)
			bg<-na.omit(bg)
            # cat('..........after NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
            
            # PR Curve
            pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
#             cat('PR_auc:',pr$auc.integral, '\n')
#             plot(pr)

            # ROC Curve    
            roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
            cat('PR_auc:',pr$auc.integral, 'ROC_auc:',roc$auc, '\n')
#             plot(roc)              
            results_by_perc[sample_count,]<-c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2) )              
        } 
	   coeff_file<-paste(datadir,"/cox_filter/method_", method,"_", 5, cancertype, missing_perc*100, ".RData",sep = "")
	   print(coeff_file)
	   save(results_by_perc, file=coeff_file)
	   
       results[perc,]<- round(apply(results_by_perc,2,mean),2)   
	   print(results)	   
       perc<-perc+1
    }
    print(results)
}



