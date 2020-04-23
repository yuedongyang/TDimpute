###compute pearson correlation for datasets imputed by different methods, Figure 3 in the paper
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
  max_pair$CpG<-as.character(max_pair$CpG)
  max_pair$gene<-as.character(max_pair$gene)

  expr_file<-paste(datadir,"/bootstrap_cancer_V1/",cancertype,"10_1.csv",sep = "") 
  print(expr_file)
  shuffle_cancer <- read.table(expr_file_1,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
  aa <-data.matrix(shuffle_cancer[,2:dim(shuffle_cancer)[2]])
  cnames<-colnames(shuffle_cancer)[2:dim(shuffle_cancer)[2]]
  rnames<- shuffle_cancer[,1]
  shuffle_cancer<-matrix(aa,nrow=dim(shuffle_cancer)[1],ncol=dim(shuffle_cancer)[2]-1,dimnames=list(rnames,cnames))
  rnames<- row.names(shuffle_cancer)
  shuffle_cancer <- shuffle_cancer[!duplicated(rnames), ] 
  
  for(missing_perc in c(0.1,0.3,0.5,0.7,0.9)){	      
	samples<-5
    coeff_matrix_1 <- matrix(nrow = samples, ncol=nrow(max_pair))
	
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
		for(i in c(1:nrow(max_pair))){
			res<-cor( x=datExpr[row.names(shuffle_cancer), max_pair[i,]$gene], y=shuffle_cancer[, max_pair[i,]$CpG] )
			max_cor[i] <- res
		}		
		max_cor[is.na(max_cor)]<-0
		coeff_matrix_1[sample_count,] <- max_cor**2
    }
	
	coeff_matrix_2 <- matrix(nrow = samples, ncol=nrow(max_pair))
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
	  for(i in c(1:nrow(max_pair))){
	    res<-cor( x=datExpr[row.names(shuffle_cancer), max_pair[i,]$gene], y=shuffle_cancer[, max_pair[i,]$CpG] )
		max_cor[i] <- res
	  }		
	  max_cor[is.na(max_cor)]<-0
	  coeff_matrix_2[sample_count,] <- max_cor**2
    }
	
	coeff_matrix_3 <- matrix(nrow = samples, ncol=nrow(max_pair))
    for(sample_count in 1:samples){
      ######################################
      expr_file<-paste(datadir,"/filled_data/TDimpute_self_",cancertype,missing_perc*100,".0_",sample_count,".csv",sep = "") 
      print(expr_file)  
	  datExpr <- read.table(expr_file,sep = ",",header=T,stringsAsFactors=FALSE,check.names=F)
	  aa <-data.matrix(datExpr[,2:dim(datExpr)[2]])
	  cnames<-colnames(datExpr)[2:dim(datExpr)[2]]
	  rnames<- datExpr[,1]
	  datExpr<-matrix(aa,nrow=dim(datExpr)[1],ncol=dim(datExpr)[2]-1,dimnames=list(rnames,cnames))
	  rnames<- row.names(datExpr)
	  datExpr <- datExpr[!duplicated(rnames), ]   
	  max_cor <- vector(mode = "numeric", length = nrow(max_pair))
	  for(i in c(1:nrow(max_pair))){
	    res<-cor( x=datExpr[row.names(shuffle_cancer), max_pair[i,]$gene], y=shuffle_cancer[, max_pair[i,]$CpG] )
		max_cor[i] <- res	
	  }		
	  max_cor[is.na(max_cor)]<-0
	  coeff_matrix_3[sample_count,] <- max_cor**2
    }

	coeff_matrix_4 <- matrix(nrow = samples, ncol=nrow(max_pair))
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
	  for(i in c(1:nrow(max_pair))){
	    res<-cor( x=datExpr[row.names(shuffle_cancer), max_pair[i,]$gene], y=shuffle_cancer[, max_pair[i,]$CpG] )
		max_cor[i] <- res
	  }
	  max_cor[is.na(max_cor)]<-0
	  coeff_matrix_4[sample_count,] <- max_cor**2
    }
	
	coeff_matrix_5 <- matrix(nrow = samples, ncol=nrow(max_pair))
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
	  for(i in c(1:nrow(max_pair))){
	    res<-cor( x=datExpr[row.names(shuffle_cancer), max_pair[i,]$gene], y=shuffle_cancer[, max_pair[i,]$CpG] )
		max_cor[i] <- res
	  }
	  max_cor[is.na(max_cor)]<-0
	  coeff_matrix_5[sample_count,] <- max_cor**2
    }
	coeff_ave_1 <- na.omit(colMeans(coeff_matrix_1))
	coeff_ave_2 <- na.omit(colMeans(coeff_matrix_2))
	coeff_ave_3 <- na.omit(colMeans(coeff_matrix_3))
	coeff_ave_4 <- na.omit(colMeans(coeff_matrix_4))
	coeff_ave_5 <- na.omit(colMeans(coeff_matrix_5))

	baseline_coeff_vec <- data.frame(value = max_pair$max_pcc, group="FULL")
	coeff_ave_1 <- data.frame(value = coeff_ave_1, group="SVD")
	coeff_ave_2 <- data.frame(value = coeff_ave_2, group="TOBMI")
	coeff_ave_3 <- data.frame(value = coeff_ave_3, group="TDimpute_self")
	coeff_ave_4 <- data.frame(value = coeff_ave_4, group="TDimpute")
	coeff_ave_5 <- data.frame(value = coeff_ave_5, group="Lasso")
	dataset <- rbind(baseline_coeff_vec, coeff_ave_1, coeff_ave_2, coeff_ave_3, coeff_ave_4, coeff_ave_5)

	coeff_file<-paste(datadir,"/coeff_matrix/corr_",samples,cancertype,missing_perc*100,"_gene_CpG_pair.RData",sep = "")
	print(coeff_file)
	save(dataset,max_pair, coeff_matrix_1,coeff_matrix_2,coeff_matrix_3,coeff_matrix_4,coeff_matrix_5, file=coeff_file)	
  }
}
