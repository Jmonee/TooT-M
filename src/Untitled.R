getPsiAndKs<- function(features)
{
  
  psi=as.integer(features/10.0001)
  K<- as.integer(((features/10.0001)-(as.integer(features/10.0001))+.01)*10)
 All<- list(psi=psi, K=K)
  return(All) 
}