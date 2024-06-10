# Ensemble Species Distribution Models

*Script is currenty being modified due to R spatial compatibility issues and dependencies on the raster package*

The SDM R script consists of 5 steps:

- Generating and formatting predictor variables
- Constraining occurrence by forest cover
- Reduce colinearity in the covariates
- Generate ESDM 
- Rescale ESDM to quantiles
- The code is for the generation of the anoa SDM, but can easily be adapted for the babirusa, any other species or for different covariates 

## 1. Preparing predictor variables
Input = raster layer (geoTiff, .tiff)
- Predictor variables are the bioclim rasters, from WorldClim version 2.1 (https://www.worldclim.org/data/worldclim21.html). There are 19 global bioclimatic rasters plus the altitude raster layer.
- Altitude is converted to topographic ruggedness index - which uncouples elevation and climate. 
- All predicator layers need to be in the same resolution and projection, and should be clipped and masked to the same geographic extent.

## 2. Constraining occurrence by forest cover
Filtered and curated occurrence data for the taxa need to be in .csv format, with one occurrence per line and a x/y coordinate column 

The occurrences need to be converted to a spatial dataframe. The points are the intersected with the forest cover vector layer (.shp), resulting in only points occurring within the forest cover layer for the generation of the model.

## 3. Testing for colinearity in the covariates
Predictor variables that are highly correlated are removed (th=3)

## 4. Generate ESDM
Using the `ensemble_modelling` function in the SSDM package, 9 algorithms are tested. The ensemble model is then generated using the algorithms that pass the 0.7 AUC threshold

The raw ESDM is then constrained by the forest cover layer - by converting any non-forested pixels to 0, giving us the habitat suitability score
HabitatSuitability = ESDMraster * ForestCoverLayer

## 5. Rescale ESDM
The ESDM is the rescaled based on the percentiles of habitat suitability scores in 20% quantile intervals which represent 5 classes (1 = low suitability, 5 = high suitability)

## Further processing - QGIS
- QGIS is then used to calculate the amounts of each habitat suitability class that can be found in the different land use categories
- This is calculated using the field calculator and overlap tool
