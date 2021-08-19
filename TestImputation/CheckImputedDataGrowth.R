#Imputation of the missing reproductive data

setwd("~/SpeciesChoice/DataPreparationSpeciesChoice/TestImputation")
SpeciesData<-readRDS(file = "~/SpeciesChoice/DataPreparationSpeciesChoice/Imputation/GrowthBeforeImpute.RData")

#load imputation package
library(mice)
library(VIM)
library(ggplot2)
library(dplyr)

source('~/SpeciesChoice/DataPreparationSpeciesChoice/DFMeansAny.R')

#prepare dataset

#identify false datapoints and remove them

SpeciesDataImp <- subset(SpeciesData, select=-species)
SpeciesDataImp$growthrateMin[SpeciesDataImp$growthrateMin<0]<-NA
SpeciesDataImp$growthrateMin[SpeciesDataImp$growthrateMin=='NaN']<-NA
SpeciesDataImp$growthrateMax[SpeciesDataImp$growthrateMax=='NaN']<-NA
SpeciesDataImp$ID <- seq.int(nrow(SpeciesDataImp))
ImputeMethods <- c('','cart','rf','norm','norm.nob','norm.boot','norm.predict','polr','polyreg','lda')

NotMissing <- filter(SpeciesDataImp,SpeciesDataImp$growthrateMin!='NA')
Loose <- NotMissing$ID


for (i in 2:3){

  #set the method for the imputation
  ImputeMethod <- ImputeMethods[i]
  for (sets in 1:2) {
    if (sets==1){
      numbersets <-50
    } else {
      numbersets<-10
    }
  MinMinimum <- vector()
  MinMean <- vector()
  MinMaximum <- vector()
  
  MaxMinimum <- vector()
  MaxMean <- vector()
  MaxMaximum <- vector()
  
  for (tries in 1:nrow(NotMissing)){
  
  SpeciesDataImp2 <- SpeciesDataImp
  SpeciesDataImp2$growthrateMin[NotMissing$ID[tries]]<- NA
  SpeciesDataImp2$growthrateMax[NotMissing$ID[tries]]<- NA
  
  imputed_Data <- mice(SpeciesDataImp2, m=numbersets, maxit = 50, method = ImputeMethod, seed = 500)
  completeData <-list()
  for (n in 1:numbersets){
    completeData <- c(completeData,list(complete(imputed_Data,n)))
  }
  
  completeData_reduced <- list()
  #Calculate range for the imputed values
  for (r in 1: length(completeData)){
    completeData_reduced[[r]] <-subset(completeData[[r]],select=growthrateMin)
  }
  
  DatasetImputedMeans <- DFMeansAny(completeData_reduced)
  DatasetImputedMeans <- data.frame(DatasetImputedMeans)

  MinMinimum[tries]=DatasetImputedMeans[,1]-DatasetImputedMeans[,2] 
  MinMean[tries]=DatasetImputedMeans[,1]
  MinMaximum[tries]=DatasetImputedMeans[,1]+DatasetImputedMeans[,2]
  
  for (r in 1: length(completeData)){
    completeData_reduced[[r]] <-subset(completeData[[r]],select=growthrateMax)
  }
  
  DatasetImputedMeans <- DFMeansAny(completeData_reduced)
  DatasetImputedMeans <- data.frame(DatasetImputedMeans)
  
  MaxMinimum[tries]=DatasetImputedMeans[,1]-DatasetImputedMeans[,2] 
  MaxMean[tries]=DatasetImputedMeans[,1]
  MaxMaximum[tries]=DatasetImputedMeans[,1]+DatasetImputedMeans[,2]
  }
  
  #Create table
  FinalTable <- data.frame(OriginalMin=NotMissing$growthrateMin,
                          ImputedMinMin=MinMinimum,
                          ImputedMinMean=MinMean,
                          ImputedMinMax=MinMaximum,
                          OriginalMax=NotMissing$growthrateMin,
                          ImputedMaxMin=MaxMinimum,
                          ImputedMaxMean=MaxMean,
                          ImputedMaxMax=MaxMaximum)
  write.csv(FinalTable,file='Growth_ImputedChecks.csv')
  write.csv(NotMissing,file='Growth_NotMissing.csv')

  }
}


#read the data back in and plot
ResultsTable<-read.csv('Growth_ImputedChecks.csv')

x<- 1:nrow(ResultsTable)
avg<- ResultsTable$ImputedMinMean
min<-ResultsTable$ImputedMinMin
max<-ResultsTable$ImputedMinMax
orig<-ResultsTable$OriginalMin

plot(x, avg,
     ylim=range(c(min, max)),
     pch=19, xlab="Measurements", ylab="Mean +/- SD",
     main="Scatter plot with std.dev error bars"
)
points(x,orig)

# hack: we draw arrows but with very special "arrowheads"
arrows(x, min, x, max, length=0.05, angle=90, code=3)

x<- 1:nrow(ResultsTable)
avg<- ResultsTable$ImputedMaxMean
min<-ResultsTable$ImputedMaxMin
max<-ResultsTable$ImputedMaxMax
orig<-ResultsTable$OriginalMax

plot(x, avg,
     ylim=range(c(min, max)),
     pch=19, xlab="Measurements", ylab="Mean +/- SD",
     main="Scatter plot with std.dev error bars"
)
points(x,orig)

# hack: we draw arrows but with very special "arrowheads"
arrows(x, min, x, max, length=0.05, angle=90, code=3)
