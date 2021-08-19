#Imputation of the missing reproductive data

setwd("~/SpeciesChoice/DataPreparation")
SpeciesData<-readRDS(file = "SkeletalBeforeImpute.RData")

#load imputation package
library(mice)
library(VIM)
library(ggplot2)
library(dplyr)

#prepare dataset

#identify false datapoints and remove them

SpeciesDataImp <- subset(SpeciesData, select=-species)
SpeciesDataImp$skeletaldensityMin[SpeciesDataImp$skeletaldensityMin<0]<-NA
SpeciesDataImp$ID <- seq.int(nrow(SpeciesDataImp))
ImputeMethods <- c('','cart','rf','norm','norm.nob','norm.boot','norm.predict','polr','polyreg','lda')

NotMissing <- filter(SpeciesDataImp,SpeciesDataImp$skeletaldensityMin!='NA')
Loose <- round(0.1*nrow(NotMissing))
RandomNumbers <- floor(runif(Loose, min=1, max=nrow(NotMissing)))
RowsLoose <- NotMissing$ID[RandomNumbers]

SpeciesDataImp2 <- subset(SpeciesDataImp, select=-ID)
LostData <- data.frame(LostRow=RowsLoose,Min=SpeciesDataImp$skeletaldensityMin[RowsLoose],Max=SpeciesDataImp$skeletaldensityMax[RowsLoose])

SpeciesDataImp2$skeletaldensityMin[RowsLoose]<-NA
SpeciesDataImp2$skeletaldensityMax[RowsLoose]<-NA

for (i in 2:3){

  #set the method for the imputation
  ImputeMethod <- ImputeMethods[i]
  numbersets <-50
  
  imputed_Data <- mice(SpeciesDataImp2, m=numbersets, maxit = 50, method = ImputeMethod, seed = 500)
  completeData <-list()
  for (n in 1:numbersets){
    completeData <- c(completeData,list(complete(imputed_Data,n)))
  }
  
  
  FileName <- paste('Skeletal_CheckImputedDatab',i,'.RData',sep="")
  save(SpeciesDataImp2,imputed_Data,completeData, file=FileName)
  
}


source('~/SpeciesChoice/DataPreparation/DFMeansAny.R')

for (i in 2:3){
  
  
  #load the appropriate method imputation result
  FileName <- paste('Skeletal_CheckImputedDatab',i,'.RData',sep="")
  load(FileName)
  
  
  #Calculate range for the imputed values
  for (r in 1: length(completeData)){
    completeData_reduced[[r]] <-subset(completeData[[r]],select=skeletaldensityMin)
  }
  
  DatasetImputedMeans <- DFMeansAny(completeData_reduced)
  DatasetImputedMeans <- data.frame(DatasetImputedMeans)

  MinMinimum=DatasetImputedMeans[,1]-DatasetImputedMeans[,2] 
  MinMean=DatasetImputedMeans[,1]
  MinMaximum=DatasetImputedMeans[,1]+DatasetImputedMeans[,2]
  
  for (r in 1: length(completeData)){
    completeData_reduced[[r]] <-subset(completeData[[r]],select=skeletaldensityMax)
  }
  
  DatasetImputedMeans <- DFMeansAny(completeData_reduced)
  DatasetImputedMeans <- data.frame(DatasetImputedMeans)
  
  MaxMinimum=DatasetImputedMeans[,1]-DatasetImputedMeans[,2] 
  MaxMean=DatasetImputedMeans[,1]
  MaxMaximum=DatasetImputedMeans[,1]+DatasetImputedMeans[,2]
  
  #Create table
  FinalTable <- data.frame(OriginalMin=LostData$Min,
                          ImputedMinMin=MinMinimum[LostData$LostRow+1],
                          ImputedMinMean=MinMean[LostData$LostRow+1],
                          ImputedMinMax=MinMaximum[LostData$LostRow+1],
                          OriginalMax=LostData$Max,
                          ImputedMaxMin=MaxMinimum[LostData$LostRow+1],
                          ImputedMaxMean=MaxMean[LostData$LostRow+1],
                          ImputedMaxMax=MaxMaximum[LostData$LostRow+1])
  

}



