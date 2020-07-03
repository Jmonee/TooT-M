##################################################
## Name: integrative.R
## Date: 25th Oct 2018
## Author: Munira Alballa
##################################################
require(seqinr)
library("Biostrings")
library("stringr")
require(protr)
library(seqinr)
library(evclass)
library(ISLR)
library(e1071)
library(caret)
library(hash)
library("stringr")
require("rlist")


pop.sd<-function(x){sqrt(sum((x-mean(x))^2)/length(x))}
standardized<-function(x,rmean,rsd){((x-rmean)/rsd)}
psiblastSeq<- function(seq, start.pos = 1L, end.pos = nchar(seq), 
                       psiblast.path = NULL,  
                       database.path = NULL, silent = TRUE, 
                       evalue = 10L, numOfIterations=3, output.path=path){
  
  
  if (Sys.which('psiblast') == '' & is.null(psiblast.path))
    stop('Please install blastp (included in the NCBI BLAST+) or specify psiblast.path.')
  
  blastp.path = if (!is.null(psiblast.path)) psiblast.path else Sys.which('psiblast')
  
  if (is.null(database.path)) stop('Must specify the database (a FASTA file) path')
  if (is.null(output.path)) stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  
  # Basic parameters for Blastp
  tmp = tempfile('Blastp2')
  queryFasta = paste0(tmp, '.fasta')
  tabularfile= paste0(tmp, '.txt')
  querySeq = Biostrings::AAStringSet(as.character(seq))
  Biostrings::writeXStringSet(querySeq, queryFasta)          
  # Additional parameters for Blastp
  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }
  # Run Blastp
  cmdblastp = paste(
    paste0(shQuote(blastp.path),
           ' -db ', shQuote(database.path),
           ' -query ', shQuote(queryFasta),  ' -num_iterations ',numOfIterations ,' -evalue ',evalue,
           
           " -out psiblastDsbAOut -out_pssm PSSMDsbA -out_ascii_pssm ascii_mtx_file"))
  
  print(cmdblastp)
  if (silent == TRUE) system(cmdblastp, ignore.stdout = TRUE) else system(cmdblastp)      
  
  
}
getTOPCONS2_predictions<- function(topconsresultsFile)
{
  fc= file(topconsresultsFile)
    mylist <-readLines(fc)
  close(fc)
  
  startOf<- grep("Sequence name: ",mylist)
  names<- sub("Sequence name: ","",mylist[startOf])
  names<-  sub("\\|.*","",sub(".+?\\|","",names))
  statistics <- matrix(ncol=100,nrow= length(names))
  
  stratfrom=1
  endfrom=1
  for(i in 1: length(names))
  {
    if(i== length(names))
    {
      endfrom= length(mylist)
    }else{
      endfrom=startOf[i+1]
      
    }
    
    seqinfo=mylist[stratfrom:endfrom]
    s<-  grep("TOPCONS predicted topology:",seqinfo)+1
    e<-  (grep("OCTOPUS predicted topology:",seqinfo, fixed = T)-2) [1]
    
    
    
    TMS<- which(unlist(strsplit(seqinfo[s:e][1],""))=="M")
    nTMS=0
    from=0
    to=0
    if(length(TMS) !=0)
    {
      ranges<- c()
      isStart=T
      for(l in 1:(length(TMS)-1))
      {
        if(TMS[l]!=(TMS[l+1]-1)){
          to=TMS[l]
          ranges<- c(ranges,to)
          isStart=T
        }else{
          if(isStart){
            nTMS=nTMS+1
            from=TMS[l]
            ranges<- c(ranges,from)
            isStart=F
          }
        }
        
        
      }
      ranges<- c(ranges, TMS[length(TMS)])
      info<- c(nTMS,ranges)
      
      statistics[i,1:length(info)]<- info
      
    }
    else{
      statistics[i,1]<- 0
    }
    stratfrom=startOf[i+1]
    
    
  }
  return(cbind.data.frame(names,ifelse(statistics[,1]>0,1,2)))
}
  
  
  
  


extractPsePSSM<- function(fastafile) #calculates Pse PSSM for lambda in (0,1,...,48,49)
{
  # Step 1
  #read the fasta file
  allData<- read.fasta(fastafile,seqtype = "AA", as.string = T)
  #extract the names of sequences
  SeqNames<- sub("\\|.*","",sub(".+?\\|","",names(allData)))
  

  for (x in c(1:length(SeqNames)) )
  {

    dirName= paste0(intermediateFiles,"/",SeqNames[x],"/")
    dir.create(dirName, showWarnings = TRUE, recursive = FALSE, mode = "0777")

    seq<- allData[x]


    setwd(dirName)
    psiblastSeq(seq, database.path=paste0(dbpath,"SwissOct18.fasta")  ,evalue=0.001, numOfIterations=3, output.path= dirName)
    # 1- reading the PSSM file
    temp= read.table(paste0(dirName,"ascii_mtx_file"), skip=3,fill=TRUE, stringsAsFactors = FALSE)
    lastRow<- which(temp$V1=="K")-1
    LastCol<- 22
    temp2=temp[1:lastRow, 1:LastCol]
    OriginalPSSM<- data.matrix(sapply(temp[1:lastRow, 3:LastCol], function(x) as.numeric(x)))
    st<-scale(OriginalPSSM[1,])

    colnames(OriginalPSSM)<- c( "A" , "R" , "N" , "D" , "C" , "Q" , "E" , "G" , "H" , "I" , "L" , "K" , "M" , "F" , "P" , "S" , "T" , "W" , "Y" , "V")
    ## these scores are then followed by a standardization procedure
    standardizedPSSM<- OriginalPSSM
    for( i in 1:length(OriginalPSSM[,1]) )#until L
    {
      #compute mean across the 20 aa
      rmean= mean(OriginalPSSM[i,])
      rsd=pop.sd(OriginalPSSM[i,])
      standardizedPSSM[i,]<- standardized(OriginalPSSM[i,],rmean,rsd)
    }



    # step 2 -> make the PSSM descriptor become a size-uniform matrix
    correaltionMatrix<- matrix(NA,nrow = 50, ncol=40)
    colnames(correaltionMatrix)<- c(rep( c( "A" , "R" , "N" , "D" , "C" , "Q" , "E" , "G" , "H" , "I" , "L" , "K" , "M" , "F" , "P" , "S" , "T" , "W" , "Y" , "V"),2))
    rownames(correaltionMatrix)<- c(paste0("PSSM",0:49))

    PSSM<- as.matrix(standardizedPSSM[,])
    length(PSSM[,1])/sum(PSSM[,1])
    #make the PSSM descriptor become a size-uniform matrix
    PSSMuniform<- apply(PSSM,2, mean, na.rm=TRUE)
    correaltionMatrix[1,1:20]<- PSSMuniform #first coorelation coefficient is zero  (Xi = 0) one 20-D (dimensional) vector (n = 0)
    for( Xi in 1:49 )
    {
      L= length(PSSM[,1])
      denominator= L- Xi

      GvectorforXi<- rep(0, 20)
      for(j in 1:20 )#length(PSSM[1,]) until w0
      {
        G.Xi.j= paste0("G.",Xi,".",j)

        sum=0
        for( i in 1:denominator)
        {
          print(((PSSM[i,j]-PSSM[i+Xi,j])^2))
          sum=sum+ ifelse(is.nan(((PSSM[i,j]-PSSM[i+Xi,j])^2)),0,(PSSM[i,j]-PSSM[i+Xi,j])^2)
        }
        GvectorforXi[j]<- sum/denominator
      }
      correaltionMatrix[Xi+1,1:20]<- PSSMuniform
      correaltionMatrix[Xi+1,21:40]<- GvectorforXi


    }

    write.csv(correaltionMatrix,paste0(dirName,"AllFeaturesPSSMuniform.csv"))
  }

  
  # create 50 Pse-PSSM (one for each Xi) in the compostion folder:
  assign(paste0("pse",1,"PSSM"), data.frame(matrix(ncol = 20, nrow =0)))
  for(i in 2:50)
  {
    assign(paste0("pse",i,"PSSM"), data.frame(matrix(ncol = 40, nrow =0)))
  }
  
  
  
  
  for( n in SeqNames)
  {
    dirName= paste0(intermediateFiles,n,"/")
    AllFeaturesPSSMuniform<- read.csv(paste0(dirName,"AllFeaturesPSSMuniform.csv"), stringsAsFactors = FALSE)[,-1]
    pse1PSSM<- rbind(pse1PSSM,unname(c(AllFeaturesPSSMuniform[i,1:20])))
    
    for(i in 2:50)
    {
      assign(paste0("pse",i,"PSSM"), rbind(get(paste0("pse",i,"PSSM")),unname(c(AllFeaturesPSSMuniform[i,]))))
    }
    
  }
  
  
  
  
  
  for(i in 1:50)
  {
    write.csv(get(paste0("pse",i,"PSSM")),file = paste0(compostions,"Testing/",paste0("pse",(i-1),"PSSM.csv")), row.names  = SeqNames)
  }
  
}
