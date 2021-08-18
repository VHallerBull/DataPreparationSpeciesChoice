setwd("~/SpeciesChoice/DataPreparationSpeciesChoice/DataDownload")

##Species list identification
#load species from Corals of the world that are present in the GBR
COTW<- read.csv("SpeciesList.csv", header =FALSE)

#load species taxonomy from Josh/Mike
TaxonomyDatabase <- read.csv("species_choice-master/data/taxonomy.csv", header =FALSE)

#match GBR speceis with the new taxonomy
NewSpecies <- data.frame(species_COTW=character(),
                         species_old=character(),species_new=character(),
                         family_old=character(),family_new=character(),
                         major_clade=character(),
                         stringsAsFactors=FALSE)


#for (i in 1:nrow(COTW)){
#  Row<-TaxonomyDatabase[TaxonomyDatabase$species_old==paste(COTW[i,1],COTW[i,2]),]
#  if (nrow(Row)==0) {
#    Row<-TaxonomyDatabase[TaxonomyDatabase$species_new==paste(COTW[i,1],COTW[i,2]),]
#    if (nrow(Row)==0) {

#      Row[1,]<-c('NA','NA','NA','NA','NA')
#    }
#  }
#  NewSpecies[i,1]<-paste(COTW[i,1],COTW[i,2])
#  NewSpecies[i,2:6]<-Row
#
#}

#save species list
#write.csv(NewSpecies,'NewSpeciesList.csv')

NewSpecies<- read.csv("NewSpeciesList.csv", header =TRUE)

#for now subset to exclude species not found in Josh/Mike taxonomy
MissingSpecies <-NewSpecies[NewSpecies$species_old=="NA",]
NewSpecies<-NewSpecies[is.na(NewSpecies$species_old)==FALSE,-1]


##Go through each trait and identify data from coraltraits.org
SpeciesData <- data.frame(species=character(),family=character())
SpeciesData[c(1:nrow(NewSpecies)),1] <- NewSpecies$species_new
SpeciesData[c(1:nrow(NewSpecies)),2] <- NewSpecies$family_new

SpeciesData_old <- data.frame(species=character(),family=character(),
                              skeletaldensityMin=double(),skeletaldensityMax=double())
SpeciesData_old[c(1:nrow(NewSpecies)),1] <- NewSpecies$species_old
SpeciesData_old[c(1:nrow(NewSpecies)),2] <- NewSpecies$family_new


#function that extracts quantitative traits as maximum and minimum from dataframe
DownloadTrait <- function(data,Species){
  TraitData <- data.frame(Min=double(), Max=double()) 
  for (i in 1:nrow(Species)){
    Trait <- data[data$specie_name==Species[i,1],]
    
    if (nrow(Trait)==0){
      Min <- NA
      Max <- NA
    } else {
      Traits_raw<-Trait[Trait$value_type=="raw_value",]
      Traits_other<-Trait[Trait$value_type!="raw_value",]
      Traits_other<-Traits_other[is.na(Traits_other$value_type)==FALSE,]
      if (nrow(Traits_raw)==0){
        Min_calc<-NA
        Max_calc<-NA
      } else {
        Min_calc<-min(Traits_raw$value[is.na(Traits_raw$value)==FALSE])
        Max_calc<-max(Traits_raw$value[is.na(Traits_raw$value)==FALSE])
      }
      if (nrow(Traits_other)>0) {
        for (m in 1:nrow(Traits_other)){
          Trait_ind<-Trait[m,]
          if (Trait_ind$value_type=='mean'){
            Mean <- as.numeric(Trait_ind$value[is.na(Trait_ind$value)==FALSE])
            if (is.na(Trait_ind$precision_type)==FALSE){
              if (Trait_ind$precision_type=="standard_deviation"){
                SD <- Trait_ind$precision[is.na(Trait_ind$precision)==FALSE]
                Min_calc[m+1] <- Mean - (3*SD)
                Max_calc[m+1] <- Mean + (3*SD)
              } else if (Trait_ind$precision_type=="standard_error") {
                SD <- Trait_ind$precision * sqrt(Trait_ind$replicates)
                Min_calc[m+1] <- Mean - (3*SD)
                Max_calc[m+1] <- Mean + (3*SD)
              }
            } else {
              Min_calc[m+1] <- Mean
              Max_calc[m+1] <- Mean
            }
          } else if (Trait_ind$value_type=="median"){
            Precision <- Trait_ind$precision[is.na(Trait_ind$precision)==FALSE]
            if (length(Precision)>0){
              Min_calc[m+1] <- Trait_ind$value - 0.5*Precision
              Max_calc[m+1] <- Trait_ind$value - 0.5*Precision
            } else {
              Min_calc[m+1] <- Trait_ind$value 
              Max_calc[m+1] <- Trait_ind$value
            }
          } 
        }}
      Min <- mean(as.numeric(Min_calc[is.na(Min_calc)==FALSE]))
      Max <- mean(as.numeric(Max_calc[is.na(Max_calc)==FALSE]))
    }
    
    TraitData[i,1] <- Min
    TraitData[i,2] <- Max
  }
  return(TraitData)
}  

#function that extracts qualitative traits  from dataframe
DownloadTraitQual <- function(data,Species){
  TraitData <- data.frame(Cat=character()) 
  for (i in 1:nrow(Species)){
    Trait <- data[data$specie_name==Species[i,1],]
    if (nrow(Trait)>0){
      TraitData[i,1] <- Trait$value[1]
    } else {
      TraitData[i,1] <- NA
    }
  }
  return(TraitData)
}

#skeletal density

data <- read.csv("https://coraltraits.org/traits/61.csv", as.is=TRUE)
data <- data[data$trait_id==61,]

TraitData_new <- DownloadTrait(data,SpeciesData)
TraitData_old <- DownloadTrait(data,SpeciesData_old)
TraitData <- data.frame(Min=double(), Max=double()) 
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,c(1,2)]<-TraitData_old[i,c(1,2)]
  } else {TraitData[i,c(1,2)]<-TraitData_new[i,c(1,2)]
  }}

SpeciesData$skeletaldensityMin <- TraitData$Min
SpeciesData$skeletaldensityMax <- TraitData$Max


#larval swimming speed

data <- read.csv("https://coraltraits.org/traits/169.csv", as.is=TRUE)
data <- data[data$trait_id==169,]

TraitData_new <- DownloadTrait(data,SpeciesData)
TraitData_old <- DownloadTrait(data,SpeciesData_old)
TraitData <- data.frame(Min=double(), Max=double()) 
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,c(1,2)]<-TraitData_old[i,c(1,2)]
  } else {TraitData[i,c(1,2)]<-TraitData_new[i,c(1,2)]
  }}

SpeciesData$larvalswimmingMin <- TraitData$Min
SpeciesData$larvalswimmingMax <- TraitData$Max


#colony shape factor

data <- read.csv("https://coraltraits.org/traits/137.csv", as.is=TRUE)
data <- data[data$trait_id==137,]


TraitData_new <- DownloadTraitQual(data,SpeciesData)
TraitData_old <- DownloadTraitQual(data,SpeciesData_old)
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,1]<-TraitData_old[i,1]
  } else {TraitData[i,1]<-TraitData_new[i,1]
  }}
SpeciesData$colonyshape <- TraitData$Min


#substrate attachement
data <- read.csv("https://coraltraits.org/traits/103.csv", as.is=TRUE)
data <- data[data$trait_id==103,]


TraitData_new <- DownloadTraitQual(data,SpeciesData)
TraitData_old <- DownloadTraitQual(data,SpeciesData_old)
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,1]<-TraitData_old[i,1]
  } else {TraitData[i,1]<-TraitData_new[i,1]
  }}
SpeciesData$substrateattachment <- TraitData$Min

#wave exposure
data <- read.csv("https://coraltraits.org/traits/144.csv", as.is=TRUE)
data <- data[data$trait_id==144,]


TraitData_new <- DownloadTraitQual(data,SpeciesData)
TraitData_old <- DownloadTraitQual(data,SpeciesData_old)
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,1]<-TraitData_old[i,1]
  } else {TraitData[i,1]<-TraitData_new[i,1]
  }}
SpeciesData$waveexposure <- TraitData$Min


#wave exposure preference
data <- read.csv("https://coraltraits.org/traits/96.csv", as.is=TRUE)
data <- data[data$trait_id==96,]


TraitData_new <- DownloadTraitQual(data,SpeciesData)
TraitData_old <- DownloadTraitQual(data,SpeciesData_old)
for (i in 1:nrow(TraitData_new)){
  if (is.na(TraitData_new[i,1])) {
    TraitData[i,1]<-TraitData_old[i,1]
  } else {TraitData[i,1]<-TraitData_new[i,1]
  }}
SpeciesData$waveexposurepreference <- TraitData$Min


#write.table(SpeciesData,"ReproductiveBeforeImpute.csv")
saveRDS(SpeciesData, file = "SkeletalBeforeImpute.RData")
