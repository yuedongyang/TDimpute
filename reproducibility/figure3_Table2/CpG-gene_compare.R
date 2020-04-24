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


###############plot summary of 16 cancers
###############top 100 CpG-gene pairs
options(stringsAsFactors=F)
samples<-5
dataset_all<-data.frame()
for(cancertype in c('LUSC', 'KIRC', 'CESC', 'STAD', 'SARC', 'COAD','KIRP', 'LUAD', 'BLCA', 'BRCA','HNSC','LGG_','PRAD','THCA','SKCM', 'LIHC')){
  coeff_file<-paste("/data/coeff_matrix/coeff_matrix_base",cancertype,"_gene_CPG_pair.RData",sep = "")
  print(coeff_file)
  load(file=coeff_file)
  max_pair$gene<-as.character(max_pair$gene)
  #sig_gene <- max_pair[max_pair$max_pcc**2>=0.5, 'gene']
  sig_gene <- max_pair[order(max_pair$max_pcc**2, decreasing = TRUE), ]
  sig_gene <- sig_gene[1:100, 'gene']  
  cat('sig_gene length:',length(sig_gene),'\n')
  
  dataset_by_cancer<-data.frame()
  for(miss_rate in c(0.1,0.3,0.5,0.7,0.9)){	 
    coeff_file<-paste("/data/coeff_matrix/corr_",samples,cancertype,miss_rate*100,"_gene_CPG_pair.RData",sep = "")
    print(coeff_file)
    load(file=coeff_file)
    method_0 <- mean(max_pair[max_pair$gene %in% sig_gene, 'max_pcc'])
    method_0 <- data.frame(correlation=method_0, method='FULL', sample_count=0)
    
    method_tmp <- data.frame(coeff_matrix_1,stringsAsFactors = FALSE)
    colnames(method_tmp)<-as.character(max_pair$gene)
    method_1<-rowMeans(method_tmp[,colnames(method_tmp) %in% sig_gene])
    method_1 <- data.frame(correlation=method_1, method='SVD', sample_count=seq(5)) 
    
    method_tmp <- data.frame(coeff_matrix_2,stringsAsFactors = FALSE)
    colnames(method_tmp)<-as.character(max_pair$gene)
    method_2<-rowMeans(method_tmp[,colnames(method_tmp) %in% sig_gene])
    method_2 <- data.frame(correlation=method_2, method='TOBMI', sample_count=seq(5))
    
    method_tmp <- data.frame(coeff_matrix_3,stringsAsFactors = FALSE)
    colnames(method_tmp)<-as.character(max_pair$gene)
    method_3<-rowMeans(method_tmp[,colnames(method_tmp) %in% sig_gene])
    method_3 <- data.frame(correlation=method_3, method='TDimpute-self', sample_count=seq(5))
    
    method_tmp <- data.frame(coeff_matrix_4,stringsAsFactors = FALSE)
    colnames(method_tmp)<-as.character(max_pair$gene)
    method_4<-rowMeans(method_tmp[,colnames(method_tmp) %in% sig_gene])      
    method_4 <- data.frame(correlation=method_4, method='TDimpute', sample_count=seq(5))
    
    method_tmp <- data.frame(coeff_matrix_5,stringsAsFactors = FALSE)
    colnames(method_tmp)<-as.character(max_pair$gene)
    method_5 <- rowMeans(method_tmp[,colnames(method_tmp) %in% sig_gene])      
    method_5 <- data.frame(correlation=method_5, method='Lasso', sample_count=seq(5))    
    
    dataset<-rbind(method_0,method_1,method_2,method_5, method_3,method_4)
    dataset$missing <- paste(miss_rate*100,'%',sep = "")
    dataset$cancertype <- cancertype
    
    if(miss_rate==0.1){ #only keep Full data in the 10% label 
      dataset[dataset$method=='FULL','missing'] <- "FULL"
    }else{
      dataset <- dataset[dataset$method!='FULL',]  ##remove the Full data in other labels (i.e., 30%-90%)
    }
    dataset_by_cancer <- rbind(dataset_by_cancer,dataset)
  }
  dataset_all <- rbind(dataset_all, dataset_by_cancer)
}

dataset_cancers<-data.frame()
dataset_mean_cancer_a0<-data.frame()
for(miss_rate in c('10%','30%','50%','70%','90%')){	  
  for(method in c('SVD','TOBMI','Lasso','TDimpute-self','TDimpute')){
    for(sample_count in 1:samples){
      ##mean of 16 cancers
      dataset<-data.frame(method=method, sample_count=sample_count,missing=miss_rate)
      aa<-mean(dataset_all[dataset_all$sample_count==sample_count & dataset_all$method==method & dataset_all$missing==miss_rate,'correlation'] )
      dataset$correlation<-aa
      dataset_cancers <- rbind(dataset_cancers, dataset)
    }
    ##mean by sample_count
    dataset<-data.frame(method=method, missing=miss_rate)
    dataset$correlation <- mean(dataset_cancers[dataset_cancers$method==method&dataset_cancers$missing==miss_rate,'correlation'])
    dataset_mean_cancer_a0 <- rbind(dataset_mean_cancer_a0, dataset)
  }
}

average_all<-data.frame(matrix(nrow = 5, ncol=5))
colnames(average_all)<-c('SVD','TOBMI','Lasso','TDimpute-self','TDimpute')
row.names(average_all)<-c('10%','30%','50%','70%','90%')
j<-0
for(method in c('SVD','TOBMI','Lasso','TDimpute-self','TDimpute')){
  j<-j+1
  i<-0
  for(miss_rate in c('10%','30%','50%','70%','90%')){	
    i<-i+1
    average_all[i,j] <- dataset_mean_cancer_a0[dataset_mean_cancer_a0$method==method&dataset_mean_cancer_a0$missing==miss_rate, 'correlation']
  }
}
print(average_all)