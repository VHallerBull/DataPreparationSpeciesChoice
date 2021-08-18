---
output: html_document
title: "ReadMe"
author: "Vanessa Haller-Bull"
---


# DataPreparationSpeciesChoice
collection and computations to prepare data for the species choice analysis


# 1. Step - DataDownload

The data for the traits of concern is downloaded from www.coraltraits.org. For quantitative traits we download the average +/- standard deviation. 

The main code to complete the data download is "DataPrep.R", it utilises the inputs on species of concern as well as their taxonomy provided by Josh.

The output is saves as "SpeciesData.csv" and contains the following traits:

Table 1: Number of missing values for the data downloaded 
Trait | Missing Data
------|-------------
skeletaldensityMin | 382
skeletaldensityMax | 382
growthrateMin | 357
growthrateMax | 357
corallitewidthMin | 50
corallitewidthMax | 20
colonydiameterMin | 254
colonydiameterMax | 254
polypfecundityMin | 403
polypfecundityMax | 403
eggsizeMin | 411
eggsizeMax | 411
photosynthesisMin | 408
photosynthesisMax | 408
BRI | 202

<<<<<<< HEAD
```{r table1, echo=FALSE, results='asis'}
library(knitr)
Trait<-c('skeletaldensityMin','skeletaldensityMax','growthrateMin','growthrateMax','corallitewidthMin','corallitewidthMax','colonydiameterMin','colonydiameterMax','polypfecundityMin','polypfecundityMax','eggsizeMin','eggsizeMax','photosynthesisMin','photosynthesisMax','BRI')
MissingData<-c(382,382,357,357,50,20,254,254,403,403,411,411,408,408,202)
Table<-data.frame(Trait,MissingData)
kable(Table, caption='Table 1: Number of missing values for the data downloaded')
```

	
=======
>>>>>>> f653b7fb5ca2cc2da26494af2076f95ab0029458
