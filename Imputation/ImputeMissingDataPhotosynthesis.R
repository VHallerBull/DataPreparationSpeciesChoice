#Imputation of the missing reproductive data

setwd("~/SpeciesChoice/DataPreparationSpeciesChoice/Imputation")
SpeciesData<-readRDS(file = "~/SpeciesChoice/DataPreparationSpeciesChoice/DataDownload/PhotosynthesisBeforeImpute.RData")

#load imputation package
library(mice)
library(VIM)
library(ggplot2)

#prepare dataset
SpeciesDataImp <- subset(SpeciesData, select=-species)
ImputeMethods <- c('','cart','rf','norm','norm.nob','norm.boot','norm.predict','polr','polyreg','lda')
png(file="Photosynthesis_Patternplot.png",width=1000)
mice_plot <- aggr(SpeciesDataImp, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(SpeciesData), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))
dev.off

for (i in 1:3){
  if (i > 3) {
    #prepare dataset
    SpeciesDataImp <- subset(SpeciesData, select=c(-species,-family,-sexualsystem,-symbiodinium))
  }
  
  #set the method for the imputation
  ImputeMethod <- ImputeMethods[i]
  numbersets <-10

  imputed_Data <- mice(SpeciesDataImp, m=numbersets, maxit = 50, method = ImputeMethod, seed = 500)
  completeData <-list()
  for (n in 1:numbersets){
    completeData <- c(completeData,list(complete(imputed_Data,n)))
  }

 
  FileName <- paste('Photosynthesis_ImputedData',i,'.RData',sep="")
  save(SpeciesDataImp,imputed_Data,completeData, file=FileName)

}

source('~/SpeciesChoice/DataPreparationSpeciesChoice/DFMeansAny.R')


for (i in 1:3){
  
  #load the appropriate method imputation result
  FileName <- paste('Photosynthesis_ImputedData',i,'.RData',sep="")
  load(FileName)
  
  
  for (n in 1:numbersets){
    DataName <- paste('completeData',n,sep="")
    
    #Plot a graph of the missing data for each dataset of the method
    png(file=paste(DataName,"_Pattern.png"))
    mice_plot_completed <- aggr(completeData[[n]], col=c('navyblue','yellow'),
                              numbers=TRUE, sortVars=TRUE,
                              labels=names(SpeciesData), cex.axis=.7,
                              gap=3, ylab=c("Missing data","Pattern"))
    dev.off()
    
    #Extract the variables that have been completely imputed
    VariablesImputed <- mice_plot_completed$missings$Variable[mice_plot_completed$missings$Count==0]
    VariablesImputed <- VariablesImputed[-1]
    
    
  }
  if (length(VariablesImputed)>0){
    completeData_reduced <- list()
    for (r in 1: length(completeData)){
      completeData_reduced[[r]] <-completeData[[r]][,VariablesImputed]
    }
         
  DatasetImputedMeans <- DFMeansAny(completeData_reduced)
  DatasetImputedMeans <- data.frame(DatasetImputedMeans)
  
  for (p in 1:length(VariablesImputed)) {
    png(file=paste(FileName,"_",VariablesImputed[p],"_plot.jpg",sep=""),width=1000)
    print(ggplot(DatasetImputedMeans,aes(c(1:nrow(DatasetImputedMeans)),DatasetImputedMeans[,p]))+
            geom_point(shape = 18, color = "#FC4E07", size = 2)+
            geom_errorbar(aes(ymin=DatasetImputedMeans[,p]-DatasetImputedMeans[,p+length(VariablesImputed)], 
                      ymax=DatasetImputedMeans[,p]+DatasetImputedMeans[,p+length(VariablesImputed)]), 
                      width=.2,position=position_dodge(.9))) 
    dev.off()
    print(p)
    
  }
  
  }
}

 
  
 
 
 