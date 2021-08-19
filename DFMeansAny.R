DFMeansAny <- function(ListDF){
  DF <- matrix(nrow=nrow(ListDF[[1]])+1,ncol=2*ncol(ListDF[[1]]))
  for (n in 1:ncol(ListDF[[1]])){
    for (i in 1:nrow(ListDF[[1]])){
      Values <- vector()
      for (l in 1:length(ListDF)){
        Values <- c(Values,ListDF[[l]][i,n])
      }
      DF[i+1,n] <- mean(Values)
      DF[i+1,n+ncol(ListDF[[1]])] <- sd(Values)
    }
    DF[1,]<-colMeans(DF[c(2:nrow(ListDF[[1]])),],)
  }

  
  return(DF)
}