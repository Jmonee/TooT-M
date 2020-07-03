
###################################################
## Name: TooT_M
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted location of each protein sequence, membrane=1, nonmembrane=2
## Author: Munira Alballa
##################################################

getlambdaAndKs<- function(features)
{
  
  lambda=as.integer(features/10.0001)
  K<- as.integer(((features/10.0001)-(as.integer(features/10.0001))+.01)*10)
  All<- list(lambda=lambda, K=K)
  return(All) 
}

TooT_M<- function(test_fasta,mRMR_file, topconspredictions)
{
  seqs<- readFASTA(test_fasta)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  results=rep(0,length( names(seqs)))
  
  mrmrInitial<-read.csv(mRMR_file)[,-1]
  testvoting=data.frame(uniprotID=names(seqs))
  
  All<- getlambdaAndKs(mrmrInitial)
  for(lambda in unique(All$lambda))
  {
    filename= paste0("pse",lambda,"PSSM.csv")
    
    taining = read.csv(paste0(compostions,"Training/final", filename), header = TRUE,sep=",", stringsAsFactors = FALSE)[,-1]
    Trdata=taining[,c(-1)]
    Trclass=taining[,1]

    test.data= read.csv(paste0(compostions,"Testing/", filename), header = TRUE,sep=",", stringsAsFactors = FALSE)[,-1]
    param0<- EkNNinit(Trdata,Trclass)
    
    for (i in All$K[which(All$psi==lambda)])
    {
    prediction<- EkNNval(Trdata,Trclass, test.data,K=i,param=param0)
    testvoting=cbind(testvoting, prediction$ypred)
    }
  }
    
    
    for(i in 1:length( names(seqs)))
    {
      temp<-unlist(as.vector(c(unname(testvoting[i,]))))
      if(topconspredictions[i,2]==1)
      {
        temp=c(temp, rep(1,19))
      }else{
        temp=c(temp, 2)
      }
      
      results[i]<- as.numeric(names(sort(table(temp),decreasing=TRUE))[1])
      
    }
    
  return(results)
}



args <- commandArgs(trailingOnly=TRUE)

terminate <- FALSE

out <- "."
TooTMdir <- "."
db<-"./db/"
topcons2<-"."
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           query <- arg[2]
         },
         "-out"={
           out <- arg[2]
         },
         "-TooTM"={
           TooTMdir <- normalizePath(arg[2])
         },
         "-db"={
           db <- normalizePath(arg[2])
         },
         "-topcons2"={
           topcons2 <- normalizePath(arg[2])
         },
         "-help"={
           cat("TooTM v1.0 (Oct. 2019)\n")
           cat("\n")
           cat("Usage: TooTM -query=<input> [-TooTM=<TooTMdir>] [-out=<outdir>] [-db=<database path>] [-topcons2=<TOPCONS2 results path>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '",out,"'\n")
           cat("\t<TooTMdir> is the directory where the base TooT-M files are located")
           cat("\t\t<TooTMdir> defaults to '",TooTMdir,"'\n")
           cat("\t <database path> is the relative path to the database\n")
           cat("\t\t<database path> defaults to",TooTMdir,"/db/\n")
           cat("\t <TOPCONS2 results path> is the query result using TOPCONS2 (.txt format) \n")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  
  if(!exists("query")) {
    stop("-query has not been passed")
  }
  
  
  
  test_fasta <- normalizePath(path.expand(query))
  resultspath <- paste0(normalizePath(path.expand(out)),"/")
  

  
  require(seqinr)
  library("Biostrings")
  library("stringr")
  require(protr)
  library(ISLR)
  library(e1071)
  library(caret)
  library(R.utils)
  wd=normalizePath(path.expand(".")) # change the the tool directory
  
  if (isAbsolutePath(db)){
    dbpath <- db
  }else{
    dbpath <- paste0(TooTMdir, db)
  }
  

  
  
  # test_fasta<- "~/TooT-M/test.fasta"
  # resultspath<- "~/TooT-M/TooTMoutput/"
  # dbpath= "~/TooT-SC.v1.svm.11/db/"
  # TooTMdir="~/TooT-M/"
  #topcons2="~/TooT-M/TOPCONS2testresults.txt"
  
  intermediateFiles=paste0(TooTMdir,"/intermediate_files/")
  compostions=paste0(TooTMdir,"/intermediate_files/Compositions/")
  
  
  #testing data with unknown substrates
  source(paste0(TooTMdir,"/src/integrative.R"))
  
  extractPsePSSM(test_fasta)
  
  TOPCONS2Predictions<-getTOPCONS2_predictions(topcons2)

  predictions<- TooT_M(test_fasta,paste0(TooTMdir,"mRMR_candidates_index.csv"),TOPCONS2Predictions)
  
  
 
  
  # write results
  seqs<- readFASTA(test_fasta)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  print(paste0( "TooT-M output is found at: ", resultspath, "TooTMout.csv"))
  write.csv(cbind(UniProtID=names(seqs),predictions ),paste0(resultspath,"TooTMout.csv"))
  
  
}
