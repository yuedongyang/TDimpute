rm(list = ls())
print('start...')

##TOBMI function used for imputation
TOBMI <- function(x = cpg, y = exp) {
  ##Calculating the distances among un-/complete cases using auxiliary dataset
  dist.matrix <- as.matrix(dist( x ))
  
  ##Neighbors list for every uncomplete cases
  missing_num <- length(which(complete.cases(y) == F)) 
  donors <- list()
  for(i in 1:missing_num){
    donors[[i]] <- as.matrix(sort(dist.matrix[i,c(c(missing_num + 1):dim(x)[1])])[1 : floor(sqrt(dim(x)[1] - missing_num))])
  }
  
  ##Neighbors will be weighted by distance 
  donors.w<-list()
  for(i in 1:missing_num){
    donors.w[[i]]<-(1/( donors[[i]][,1]+0.000000001 ))/sum((1/( donors[[i]][,1]+0.000000001 )))
    #     donors.w[[i]]<-(1/( donors[[i]][,1] ))/sum((1/( donors[[i]][,1] )))
  }
  
  ##Imputation process
  for(j in 1:missing_num){
    as.data.frame(donors.w[[j]])->donors.pool
    row.names(donors.pool)->donors.pool$id
    y$id <- row.names(y)
    merge(donors.pool,y,by='id')->donors.candidates
    donors.candidates[,2] * donors.candidates[,3:dim(donors.candidates)[2]]->donors.calculate
    y[j,-dim(y)[2]]<-apply(donors.calculate, MARGIN = 2,sum)
  }
  imputed.data<-y[,-dim(y)[2]]
}

reduced_names = c('WT')
datadir <- '/data' 

cancer_c <- 1
for(cancertype in reduced_names){
  rmse_arr=matrix(nrow = 5, ncol=5) 
  perc <- 1
  for(missing_perc in c(0.5)){
    for(sample_count in c(1,2,3,4,5)){
	  expr_file<-paste("/data",cancertype,missing_perc*100,".0_",sample_count,"_quantiles_RSEM.csv",sep = "") 
      RNA_size = 18155 
      DNA_size = 23003
      print(expr_file)
      RDNA <- read.table(expr_file, sep = ",", header = TRUE, stringsAsFactors=FALSE, check.names=F)  
      #     row.names(RDNA)<-RDNA[,1]
      #     RDNA <- RDNA[,2:dim(RDNA)[2]] # to avoid duplicated row.names error !!! duplicate 'row.names' are not allowed
      aa <-data.matrix(RDNA[,2:dim(RDNA)[2]])
      cnames<-colnames(RDNA)[2:dim(RDNA)[2]]
      rnames<- RDNA[,1]
      RDNA<-matrix(aa,nrow=dim(RDNA)[1],ncol=dim(RDNA)[2]-1,dimnames=list(rnames,cnames))
      RDNA <- RDNA[!duplicated(rownames(RDNA)), ]
	  
      print(floor(dim(RDNA)[1]*missing_perc)+1)
      test_size <- floor(dim(RDNA)[1]*missing_perc)+1
      
      RDNA_corrupted <- RDNA
      RDNA_corrupted[1:test_size,1:RNA_size] <- NA
      
      cat('dim(RDNA_corrupted)[1]: ',dim(RDNA_corrupted)[1],'\n')
      
      exp <- as.data.frame(RDNA_corrupted[,1:RNA_size])
      cpg <- as.data.frame(RDNA_corrupted[,(RNA_size+1):dim(RDNA)[2]])
      
      ##TOBMI function will return a complete dataset imputed by TOBMI
      imputed<-TOBMI(x = cpg, y = exp)
      write.csv(imputed, file=paste(datadir,'/filled_data/TOBMI_',cancertype,missing_perc*100,".0_",sample_count,'.csv',sep = ""), quote=F)
      
      rmse = (mean(as.vector(as.matrix((imputed[1:test_size,] - RDNA[1:test_size,])^2))))^0.5
      cat('missing perc: ',missing_perc ,' rmse: ',rmse, '\n')
      rmse_arr[perc, sample_count] <- rmse            
    }
    perc <- perc + 1
  }
  cancer_c <- cancer_c + 1
  print( round(apply(rmse_arr,1,mean), 3) )
  print( round(rmse_arr,3) )
}
