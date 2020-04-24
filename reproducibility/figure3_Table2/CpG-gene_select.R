##compare the methylation-driving gene lists from different imputed datasets, Table 2 in the paper
rm(list = ls())
print('start...')
library('psych')
options(stringsAsFactors=F)

datadir <- '/data' 
reduced_names = c('BRCA', "SARC", "LUSC", "BLCA", "KIRC", "LGG_", "PRAD", "LUAD", "LIHC", "SKCM",
                 "HNSC", "CESC", "COAD", "KIRP", "THCA", "STAD")
				 
for(cancertype in reduced_names){
  coeff_file<-paste(datadir,"/coeff_matrix/coeff_matrix_base",cancertype,"_gene_CPG_pair.RData",sep = "")
  print(coeff_file)
  load(file=coeff_file)
  True_pos <- as.character( max_pair[max_pair$max_pcc>0.5 & max_pair$adj_pvalue<0.05, 'gene'] )
  True_neg <- as.character( max_pair[!(max_pair$max_pcc>0.5 & max_pair$adj_pvalue<0.05), 'gene'] )

  sort_gene <- max_pair[order(max_pair$max_pcc, decreasing = TRUE), ]
  top_n<-100
  top_gene_baseline = sort_gene[1:top_n, ]

  print(cancertype)
  expr_file_1<-paste(datadir,"/bootstrap_cancer_V1/",cancertype,"10_1.csv",sep = "") 
  print(expr_file_1)
  shuffle_cancer <- read.table(expr_file_1,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
  aa <-data.matrix(shuffle_cancer[,2:dim(shuffle_cancer)[2]])
  cnames<-colnames(shuffle_cancer)[2:dim(shuffle_cancer)[2]]
  rnames<- shuffle_cancer[,1]
  shuffle_cancer<-matrix(aa,nrow=dim(shuffle_cancer)[1],ncol=dim(shuffle_cancer)[2]-1,dimnames=list(rnames,cnames))
  rnames<- row.names(shuffle_cancer)
  shuffle_cancer <- shuffle_cancer[!duplicated(rnames), ] 
  
  for(missing_perc in c(0.1,0.3,0.5,0.7,0.9)){	      
	samples<-5
	sig_matrix_1 <- matrix(nrow = samples, ncol=6)

    for(sample_count in 1:samples){
        ######################################
        expr_file<-paste(datadir,"/filled_data/SVD_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "")
		print(expr_file)
		datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
		aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
		cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
		rnames<- datExpr[,1]
		datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
		rnames<- row.names(datExpr)
		datExpr <- datExpr[!duplicated(rnames), ]
		
		max_cor <- vector(mode = "numeric", length = nrow(max_pair))
		max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
		for(i in c(1:nrow(max_pair))){
			res<-cor.test( x=datExpr[row.names(shuffle_cancer), as.character(max_pair[i,]$gene)], y=shuffle_cancer[, as.character(max_pair[i,]$CpG)] )
			max_cor[i] <- res$estimate**2
			max_pvalue[i] <- res$p.value
		}
		max_pair$test_max_pcc<-max_cor		
		max_pair$test_adj_pvalue<-p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue))			
		
		sort_gene <- max_pair[order(max_pair$test_max_pcc, decreasing = TRUE), ]
		top_n<-100
		top_gene_test = sort_gene[1:top_n, ]			
		comm_gene<-intersect(top_gene_baseline$gene, top_gene_test$gene)
		cat('comm_gene-> ','sample_count:', sample_count, length(comm_gene), '\n')				
			
		pred_pos <- as.character( max_pair[(max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05,'gene'] )
		pred_neg <- as.character( max_pair[!((max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05),'gene'] )
		TP = length(intersect(True_pos, pred_pos))
		FP = length(intersect(True_neg, pred_pos))
		FN = length(intersect(True_pos, pred_neg))

		PPV = TP/(TP+FP)
		TPR = TP/(TP+FN)
		Fscore = (2*PPV*TPR)/(PPV+TPR)
		cat('TP, FP, FN:',TP,FP,FN, '\n')        
		cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')
			
		#### ROC curve- PR curve    
		base_cut<-max_pair
		base_cut$label <- rep(0,nrow(max_pair))
		r2_cutoff <- 0.5
		base_cut$label[base_cut$max_pcc>r2_cutoff] <- 1
		base_cut$pred_score <- max_pair$test_max_pcc

		# install.packages("PRROC")
		require(PRROC)
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
		sig_matrix_1[sample_count,] <- c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2), length(comm_gene))    
    }
	
	sig_matrix_2 <- matrix(nrow = samples, ncol=6)
    for(sample_count in 1:samples){
      ######################################
      expr_file<-paste(datadir,"/filled_data/TOBMI_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "")
      print(expr_file)
	  datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
	  aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
	  cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
	  rnames<- datExpr[,1]
	  datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
	  rnames<- row.names(datExpr)
	  datExpr <- datExpr[!duplicated(rnames), ] 

		max_cor <- vector(mode = "numeric", length = nrow(max_pair))
		max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
		for(i in c(1:nrow(max_pair))){
			res<-cor.test( x=datExpr[row.names(shuffle_cancer), as.character(max_pair[i,]$gene)], y=shuffle_cancer[, as.character(max_pair[i,]$CpG)] )
			max_cor[i] <- res$estimate**2
			max_pvalue[i] <- res$p.value
		}	
		max_pair$test_max_pcc<-max_cor		
		max_pair$test_adj_pvalue<-p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue))			
		
		sort_gene <- max_pair[order(max_pair$test_max_pcc, decreasing = TRUE), ]
		top_n<-100
		top_gene_test = sort_gene[1:top_n, ]			
		comm_gene<-intersect(top_gene_baseline$gene, top_gene_test$gene)
		cat('comm_gene-> ','sample_count:', sample_count, length(comm_gene), '\n')				
			
		pred_pos <- as.character( max_pair[(max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05,'gene'] )
		pred_neg <- as.character( max_pair[!((max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05),'gene'] )
		TP = length(intersect(True_pos, pred_pos))
		FP = length(intersect(True_neg, pred_pos))
		FN = length(intersect(True_pos, pred_neg))

		PPV = TP/(TP+FP)
		TPR = TP/(TP+FN)
		Fscore = (2*PPV*TPR)/(PPV+TPR)
		cat('TP, FP, FN:',TP,FP,FN, '\n')        
		cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')
			
		#### ROC curve- PR curve    
		base_cut<-max_pair
		base_cut$label <- rep(0,nrow(max_pair))
		r2_cutoff <- 0.5
		#base_cut$label[base_cut$r2>r2_cutoff & base_cut$p<0.05] <- 1
		base_cut$label[base_cut$max_pcc>r2_cutoff] <- 1
		#base_cut$label[base_cut$r2<=r2_cutoff] <- 0
		base_cut$pred_score <- max_pair$test_max_pcc
		#base_cut$pred_score <- 1-as.vector(coeff$p)

		# install.packages("PRROC")
		require(PRROC)
		# fg <- probs[df$label == 1]
		# bg <- probs[df$label == 0]
		fg <- base_cut$pred_score[base_cut$label == 1]
		bg <- base_cut$pred_score[base_cut$label == 0]
		cat('..........before NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
		# fg <- fg[!is.na(fg)]
		# bg <- bg[!is.na(bg)]
		#fg[is.na(fg)] <- 1 
		#bg[is.na(bg)] <- 1
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
		sig_matrix_2[sample_count,] <- c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2), length(comm_gene) ) 
    }
	
	sig_matrix_3 <- matrix(nrow = samples, ncol=6)
	#####################Autoencoder_by_cancer
    for(sample_count in 1:samples){
      ######################################
      expr_file<-paste(datadir,"/filled_data/TDimpute_self",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "")
      print(expr_file)  
	  datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
	  aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
	  cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
	  rnames<- datExpr[,1]
	  datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
	  rnames<- row.names(datExpr)
	  datExpr <- datExpr[!duplicated(rnames), ] 
	  
		max_cor <- vector(mode = "numeric", length = nrow(max_pair))
		max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
		for(i in c(1:nrow(max_pair))){
			res<-cor.test( x=datExpr[row.names(shuffle_cancer), as.character(max_pair[i,]$gene)], y=shuffle_cancer[, as.character(max_pair[i,]$CpG)] )
			max_cor[i] <- res$estimate**2
			max_pvalue[i] <- res$p.value
		}	
		max_pair$test_max_pcc<-max_cor		
		max_pair$test_adj_pvalue<-p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue))			
		
		sort_gene <- max_pair[order(max_pair$test_max_pcc, decreasing = TRUE), ]
		top_n<-100
		top_gene_test = sort_gene[1:top_n, ]			
		comm_gene<-intersect(top_gene_baseline$gene, top_gene_test$gene)
		cat('comm_gene-> ','sample_count:', sample_count, length(comm_gene), '\n')				
			
		pred_pos <- as.character( max_pair[(max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05,'gene'] )
		pred_neg <- as.character( max_pair[!((max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05),'gene'] )
		TP = length(intersect(True_pos, pred_pos))
		FP = length(intersect(True_neg, pred_pos))
		FN = length(intersect(True_pos, pred_neg))

		PPV = TP/(TP+FP)
		TPR = TP/(TP+FN)
		Fscore = (2*PPV*TPR)/(PPV+TPR)
		cat('TP, FP, FN:',TP,FP,FN, '\n')        
		cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')
			
		#### ROC curve- PR curve    
		base_cut<-max_pair
		base_cut$label <- rep(0,nrow(max_pair))
		r2_cutoff <- 0.5
		#base_cut$label[base_cut$r2>r2_cutoff & base_cut$p<0.05] <- 1
		base_cut$label[base_cut$max_pcc>r2_cutoff] <- 1
		#base_cut$label[base_cut$r2<=r2_cutoff] <- 0
		base_cut$pred_score <- max_pair$test_max_pcc
		#base_cut$pred_score <- 1-as.vector(coeff$p)

		# install.packages("PRROC")
		require(PRROC)
		# fg <- probs[df$label == 1]
		# bg <- probs[df$label == 0]
		fg <- base_cut$pred_score[base_cut$label == 1]
		bg <- base_cut$pred_score[base_cut$label == 0]
		cat('..........before NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
		# fg <- fg[!is.na(fg)]
		# bg <- bg[!is.na(bg)]
		#fg[is.na(fg)] <- 1 
		#bg[is.na(bg)] <- 1
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
		sig_matrix_3[sample_count,] <- c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2), length(comm_gene) ) 

    }

	sig_matrix_4 <- matrix(nrow = samples, ncol=6)
    for(sample_count in 1:samples){
      ######################################
      expr_file<-paste(datadir,"/filled_data/TDimpute_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "")
      print(expr_file)
	  datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
	  aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
	  cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
	  rnames<- datExpr[,1]
	  datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
	  rnames<- row.names(datExpr)
	  datExpr <- datExpr[!duplicated(rnames), ] 

		max_cor <- vector(mode = "numeric", length = nrow(max_pair))
		max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
		for(i in c(1:nrow(max_pair))){
			res<-cor.test( x=datExpr[row.names(shuffle_cancer), as.character(max_pair[i,]$gene)], y=shuffle_cancer[, as.character(max_pair[i,]$CpG)] )
			max_cor[i] <- res$estimate**2
			max_pvalue[i] <- res$p.value
		}	
		max_pair$test_max_pcc<-max_cor		
		max_pair$test_adj_pvalue<-p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue))			
		
		sort_gene <- max_pair[order(max_pair$test_max_pcc, decreasing = TRUE), ]
		top_n<-100
		top_gene_test = sort_gene[1:top_n, ]			
		comm_gene<-intersect(top_gene_baseline$gene, top_gene_test$gene)
		cat('comm_gene-> ','sample_count:', sample_count, length(comm_gene), '\n')				
			
		pred_pos <- as.character( max_pair[(max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05,'gene'] )
		pred_neg <- as.character( max_pair[!((max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05),'gene'] )
		TP = length(intersect(True_pos, pred_pos))
		FP = length(intersect(True_neg, pred_pos))
		FN = length(intersect(True_pos, pred_neg))

		PPV = TP/(TP+FP)
		TPR = TP/(TP+FN)
		Fscore = (2*PPV*TPR)/(PPV+TPR)
		cat('TP, FP, FN:',TP,FP,FN, '\n')        
		cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')
			
		#### ROC curve- PR curve    
		base_cut<-max_pair
		base_cut$label <- rep(0,nrow(max_pair))
		r2_cutoff <- 0.5
		#base_cut$label[base_cut$r2>r2_cutoff & base_cut$p<0.05] <- 1
		base_cut$label[base_cut$max_pcc>r2_cutoff] <- 1
		#base_cut$label[base_cut$r2<=r2_cutoff] <- 0
		base_cut$pred_score <- max_pair$test_max_pcc
		#base_cut$pred_score <- 1-as.vector(coeff$p)

		# install.packages("PRROC")
		require(PRROC)
		# fg <- probs[df$label == 1]
		# bg <- probs[df$label == 0]
		fg <- base_cut$pred_score[base_cut$label == 1]
		bg <- base_cut$pred_score[base_cut$label == 0]
		cat('..........before NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
		# fg <- fg[!is.na(fg)]
		# bg <- bg[!is.na(bg)]
		#fg[is.na(fg)] <- 1 
		#bg[is.na(bg)] <- 1
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
		sig_matrix_4[sample_count,] <- c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2), length(comm_gene) ) 
    }
	
	sig_matrix_5 <- matrix(nrow = samples, ncol=6)
    for(sample_count in 1:samples){
      ######################################
      expr_file<-paste(datadir,"/filled_data/Lasso_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "")
      print(expr_file)
	  datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
	  aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
	  cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
	  rnames<- datExpr[,1]
	  datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
	  rnames<- row.names(datExpr)
	  datExpr <- datExpr[!duplicated(rnames), ] 

		max_cor <- vector(mode = "numeric", length = nrow(max_pair))
		max_pvalue <- vector(mode = "numeric", length = nrow(max_pair))
		for(i in c(1:nrow(max_pair))){
			res<-cor.test( x=datExpr[row.names(shuffle_cancer), as.character(max_pair[i,]$gene)], y=shuffle_cancer[, as.character(max_pair[i,]$CpG)] )
			max_cor[i] <- res$estimate**2
			max_pvalue[i] <- res$p.value
		}	
		max_pair$test_max_pcc<-max_cor		
		max_pair$test_adj_pvalue<-p.adjust(max_pvalue, method ="fdr", n=length(max_pvalue))			
		
		sort_gene <- max_pair[order(max_pair$test_max_pcc, decreasing = TRUE), ]
		top_n<-100
		top_gene_test = sort_gene[1:top_n, ]			
		comm_gene<-intersect(top_gene_baseline$gene, top_gene_test$gene)
		cat('comm_gene-> ','sample_count:', sample_count, length(comm_gene), '\n')				
			
		pred_pos <- as.character( max_pair[(max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05,'gene'] )
		pred_neg <- as.character( max_pair[!((max_pair$test_max_pcc)>0.5 & (max_pair$test_adj_pvalue)<0.05),'gene'] )
		TP = length(intersect(True_pos, pred_pos))
		FP = length(intersect(True_neg, pred_pos))
		FN = length(intersect(True_pos, pred_neg))

		PPV = TP/(TP+FP)
		TPR = TP/(TP+FN)
		Fscore = (2*PPV*TPR)/(PPV+TPR)
		#cat('true_pos, true_neg, pred_pos, pred_neg:',length(true_pos), length(true_neg), length(pred_pos), length(pred_neg), '\n')
		cat('TP, FP, FN:',TP,FP,FN, '\n')        
		cat('PPV, TPR, Fsocre:',round(PPV,2), round(TPR,2), round(Fscore,2), '\n')
			
		#### ROC curve- PR curve    
		base_cut<-max_pair
		base_cut$label <- rep(0,nrow(max_pair))
		r2_cutoff <- 0.5
		#base_cut$label[base_cut$r2>r2_cutoff & base_cut$p<0.05] <- 1
		base_cut$label[base_cut$max_pcc>r2_cutoff] <- 1
		#base_cut$label[base_cut$r2<=r2_cutoff] <- 0
		base_cut$pred_score <- max_pair$test_max_pcc
		#base_cut$pred_score <- 1-as.vector(coeff$p)

		# install.packages("PRROC")
		require(PRROC)
		# fg <- probs[df$label == 1]
		# bg <- probs[df$label == 0]
		fg <- base_cut$pred_score[base_cut$label == 1]
		bg <- base_cut$pred_score[base_cut$label == 0]
		cat('..........before NA length:',sum(is.na(fg)),sum(is.na(bg)),'\n')
		# fg <- fg[!is.na(fg)]
		# bg <- bg[!is.na(bg)]
		#fg[is.na(fg)] <- 1 
		#bg[is.na(bg)] <- 1
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
		sig_matrix_5[sample_count,] <- c( round(PPV,2), round(TPR,2), round(Fscore,2), round(pr$auc.integral,2), round(roc$auc,2), length(comm_gene) ) 
    }
	print(colMeans(sig_matrix_1))
	print(colMeans(sig_matrix_2))
	print(colMeans(sig_matrix_3))
	print(colMeans(sig_matrix_4))
	print(colMeans(sig_matrix_5))

	coeff_file<-paste(datadir,"/coeff_matrix/PRAUC_data_",samples,cancertype,missing_perc*100,".RData",sep = "")
	print(coeff_file)
	save(sig_matrix_1, sig_matrix_2, sig_matrix_3, sig_matrix_4, sig_matrix_5, file=coeff_file)
  }
}
