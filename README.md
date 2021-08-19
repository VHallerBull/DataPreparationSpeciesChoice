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

Next, we download the extra data that will be used to impute missing datapoints for each trait. This is done through the code called "[Trait name]Data.R" and teh results are saved as "[Trait name]BeforeImputed.RData".

The following traits are used to support each of the traits in the main analysis.

Table 2: Traits used for the imputation for each variable in the final dataset
Variable in final dataset | Variables used for imputation
--------------------------|-----------------------------
Skeletal density | Larval swimming speed, Colony shape factor, Substrate attachment, Wave exposure, Wave exposure preference
Growth rate | Calcification rate, Life history strategy, Growth form, Growth form from Veron
Corallite width	 | Polyps per area, Tissue thickness, Total biomass
Colony diameter	 | Colony area, Colony shape factor, Coloniality
Reproduction	   | Polyp fecundity, Egg size, Colony fecundity, Mode of larval development, Propagule size, Sexual system, Symbiodinium
Photosynthesis	 | Symbiodinium density, Symbiodinium clade, Symbiodinium subclade, Zooxanthellate


# 2. Step - Imputation

The imputation trials different methods and repeats these multiple times to receive teh best results. The code is saved as "ImputeMissingData[Trait name].R" and the results are saved in multiple dataframes named "[Trait name]ImputedData[Numbering].RData".

# 3. Step - TestImputation

To test the imputation results, I iteratively remove one known datapoint and refit the model to this dataset and then predict this removed datapoint using the fitted model. Then I calculate the mean +/- standard deviation for each of these predicted datapoints which can then be compared to the known datapoint prior to removal.

The code to conduct this analysis is named "CheckImputedData[Trait name].R".

