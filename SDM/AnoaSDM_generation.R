###==================================================###
### Ensemble Species Distribution Modelling Pipeline ###
###==================================================###

## Code developed by: Dr Nicolas J. Deere, Durrell Institute of Conservation and Ecology (DICE), University of Kent
## Code edited and implemented by Dr Rosie Drinkwater, Ludwig Maximillians University (LMU)
## For troubleshooting, please contact: n.j.deere@kent.ac.uk and/or r.drinkwater@lmu.de

# Clear the globl environment
rm(list=ls())

# Establish working directory
if(!require(rstudioapi)){install.packages('rstudioapi')}
library( rstudioapi ) 

path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

# Read R packages required to perform operations
if(!require(pacman)){install.packages('pacman')}
pacman::p_load(flexsdm, spatialEco, terra, dplyr, readr, geodata, tidyverse)

# If you encounter an error installing the "flexsdm" and "geodata" packages, install manually using the code below
install.packages("remotes")
remotes::install_github("sjevelazco/flexsdm")
remotes::install_github("rspatial/geodata")

# Read packages
pacman::p_load(flexsdm, geodata)

## -----------------------------------------------------------------------------
## Define area of interest (.shp format)
## -----------------------------------------------------------------------------
# Read shapefile representing the distribution of the target species. Here, just the anoa
aoiA <- terra::vect("inputs/anoa_clipped/anoa_clipped.shp")

## -----------------------------------------------------------------------------
## Download/process predictor variables
## -----------------------------------------------------------------------------
# ............................
# Bioclim data
# ............................
#### ONLY NEED TO DO THIS ONCE THE FIRST TIME YOU DOWNLOAD THE BIOCLIMATIC VARIABLES #### 
# Download global bioclim variables, specifying the path where you want to store the data
# WARNING: this may take several minutes
worldclim <- worldclim_global(var = "bio", res = 2.5, path="worldclim_data/")

# Crop and mask to the area of interest - here just Sulawesi and the small islands 
bio.mosaic <- crop(worldclim, aoiA, mask = TRUE)

# Match the spatial extents
ext(bio.mosaic) <- ext(aoiA)

# Write the raster to file - so next time this can be loaded without having to download the global data again
writeRaster(bio.mosaic, "inputs/bioclim_2.5m_clipped.tiff", overwrite=TRUE)

## AFTER YOU DOWNLOAD/PROCESS/SAVE THE BIOCLIM VARIABLES - THE PROCESSED RASTERS CAN BE LOADED FOR SUBSEQUENT APPLICATIONS ##
## Substitute the code below for the processing functions on Lines 45-55 to load the bioclimatic rasters
#bio.mosaic <- rast("inputs/bioclim_2.5m_clipped.tiff")

#............................
# Elevation data
#............................
#### ONLY NEED TO DO THIS ONCE THE FIRST TIME YOU DOWNLOAD THE ELEVATION DATA #### 
# Download elevation for the whole world -  specify the path where you want to store the data
# Make sure the resolution matches the bioclim data
global_elev <- elevation_global(res=2.5, path="worldclim_data/")

# Crop and mask to the area of interest - here just Sulawesi and the small islands 
sulawesi_elev <- crop(global_elev, aoiA, mask = TRUE)

# Match the spatial extents
ext(sulawesi_elev) <- ext(aoiA)

# Write the raster to file - so next time this can be loaded without having to download the global data again
writeRaster(sulawesi_elev, "inputs/sulawesi_elevation2.5m_clipped.tiff", overwrite=TRUE)

## AFTER YOU DOWNLOAD/PROCESS/SAVE THE ELEVATION DATA - THE PROCESSED RASTERS CAN BE LOADED FOR SUBSEQUENT APPLICATIONS ##
## Substitute the code below for the processing functions on Lines 68-77 to load the elevation raster
# sulawesi_elev <- rast("inputs/sulawesi_elevation2.5m_clipped.tiff")

# Convert sulawesi_elev to topographic ruggedness (TRI) to uncouple elevation from climate
tri.mosaic <- spatialEco::tri(sulawesi_elev)

# Write the tri raster to file 
writeRaster(tri.mosaic, "inputs/trimosaic_clipped.tiff")

# Add the spatial predictor variables, to a raster stack and save (optional) 
sdm.vars_18 <- c(bio.mosaic, tri.mosaic)
# writeRaster(sdm.vars_18, "inputs/sdm.vars_18_all.tiff")

## -----------------------------------------------------------------------------
## Forest cover layer and processing 
## -----------------------------------------------------------------------------
# Read in the forest cover raster (downloaded from Global Forest Change: https://earthenginepartners.appspot.com/science-2013-global-forest)
forestcover <- rast("inputs/forest_18_worldclim_clipped.tif")

# Clip and mask the forest cover layer to our area of interest 
forestcover <- crop(forestcover, aoiA, mask=T)

# Convert the forest cover layer to a vector
# First, resample the forest cover raster to match the resolution to the environmental variables
# Use the nearest neighbour method of resampling as the values are discrete
FC_resample <- terra::resample(forestcover, sdm.vars_18, "near")

# Convert the forest cover raster to a vector of polygons
tmp <- as.polygons(FC_resample)

# To return a polygon of the forest outline, disaggregate the values and remove segmentsine
fc_poly <- disagg(tmp, segments = T)

## -----------------------------------------------------------------------------
## Subset species occurrence records to forest cover
## -----------------------------------------------------------------------------
# Read dataframe containing coordinates of all species occurrence records
anoa_unique <- read.csv("inputs/anoa_occ_UNIQUE2.csv", header=TRUE, stringsAsFactors=FALSE)

# Convert occurrence records to a spatial object
anoa_uni_sp <- vect(anoa_unique, geom=c("Long", "Lat"))

# Write the shapefile of occurrence records to file
writeVector(anoa_uni_sp, "outputs/anoa_occ", overwrite=T)

# Set the projection and extent of the occurrence record points to match the forest cover vector
crs(anoa_uni_sp) <- crs(fc_poly)
fc_ext <- ext(fc_poly)

# Subset the occurrence records to those within areas of forest
anoa_uni_sp_fc <- terra::crop(anoa_uni_sp, fc_poly)
plot(fc_poly)
points(anoa_uni_sp_fc)

## -----------------------------------------------------------------------------
## Extract the predictors to points, prepare for model
## -----------------------------------------------------------------------------
# Extract the coordinates from the filtered occurrence records
coords <- data.frame(crds(anoa_uni_sp_fc))

# Extract the predictor variables (bioclimatic variables, tri) to the occurrence record points
ex.vars <- as.data.frame(raster::extract(sdm.vars_18, coords))

# This follows the processing of the temperature variables https://worldclim.org/data/v1.4/formats.html 
# Write function to process temperature variables following guidance on the bioclim database
div.10 <- function(x, na.rm=FALSE) (x/10)

# Apply function to temperature variables (i.e. bio1, bio5, bio6 and bio7)
ex.vars2 <- ex.vars %>% mutate_at(c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", "wc2.1_2.5m_bio_7"), div.10, na.rm=TRUE)

# Add the coordinates to the extracted and processed covariates
anoa.vars <- cbind.data.frame(coords, ex.vars2)

# Rename the coordinate columns for modelling
colnames(anoa.vars)[1] <- "Long"
colnames(anoa.vars)[2] <- "Lat"

# Remove NAs that have crept into the data during extraction
anoa.vars <- na.omit(anoa.vars)

## -----------------------------------------------------------------------------
## Generate ESDM model
## -----------------------------------------------------------------------------
# ............................
# Remove colinear predictors 
# ............................
# Define and apply theta threshold - removes predictors with VIF above theta
vif_var18 <- flexsdm::correct_colinvar(sdm.vars_18, method = c("vif", th = "4"))

# Check twhich variables were removed 
vif_var18$removed_variables

# Access the table of variance inflation factors
vif_var18$vif_table

# Extract rasters of retained variables
(sdm.vars18_VIFreduced <- vif_var18$env_layer)

# Write the uncorrelated raster stack to file 
writeRaster(sdm.vars18_VIFreduced, "outputs/sdm_vars18_reduced_VIF5.tif", overwrite = TRUE)

#............................
# Ensemble distribution model
#............................

## Reformat as raster stack (from SpatRaster) 
sdm.vars18_VIFreduced_stack <- raster::stack(sdm.vars18_VIFreduced)

# Run the ensemble_modelling functions - requires as input variables:
# - SDM algorithms to implement  (here, 8 algorithms: 'GAM', 'GLM', 'SVM', 'MARS', 'GBM', 'CTA', 'ANN', 'RF'), 
# - Dataframe of occurrences records, including extracted predictors (anoa.vars) 
# - A raster stack containing the uncorrelated set of spatially-explicit explicit predictors (sdm.vars18_VIFreduced_stack)
# - An ensemble.metric - here we use AUC (area under the curve) to select which of the best SDMs to be included in the ensemble 
# - The ensemble.threshold by which to merge those best SDMs together in the ensemble
ESDM <- SSDM::ensemble_modelling(c('GAM', 'GLM', 'SVM', 'MARS', 'GBM', 'CTA', 'ANN', 'RF'), 
                                 anoa.vars, 
                                 sdm.vars18_VIFreduced_stack, rep = 10, Xcol = 'Long', Ycol = 'Lat', 
                                 ensemble.thresh = c(0.70), ensemble.metric="AUC", weight=TRUE, 
                                 verbose = FALSE)

## Extract the projected ESDM
anoa_ESDM_raw <- rast(ESDM@projection)
plot(anoa_ESDM_raw)

# Write the projected ESDM to file
writeRaster(anoa_ESDM_raw,  "outputs/anoa_RAW_ESDM.tif", overwrite=T)

#..........................
# Model evaluation
#..........................
# Assess model accuracy and save output to csv
# evaluating the full model
(ESDM.acc <- as.data.frame(ESDM@evaluation))
#write.csv(ESDM.acc, file="outputs/RESULTS_EDSM_acc.csv")

# Assess variable importance and save output to csv
# variable importance is assessed using the pearsons correlation (r) between the full model and the one without the variable 
(ESDM.var <- as.data.frame(ESDM@variable.importance))
#write.csv(ESDM.var, file="outputs/RESULTS_ESDM_varImportance.csv")

# Assess underlying SDM algorithms and save output to csv
# this provides the individual model evaluations
(ESDM.alg <- as.data.frame(ESDM@algorithm.evaluation))
#write.csv(ESDM.alg, file="outputs/RESULT_ESDM_algorithmEvaluation.csv")

## -----------------------------------------------------------------------------
# Processing the ESDM 
## -----------------------------------------------------------------------------
#...........................................
# Subset ESDM to areas of forest cover only
#...........................................
# Mask the ESDM to forest cover (extents need to match)
# This can be achieved by multiplying the ESDM by the *resampled* forest cover raster
FC_habitatSui <- anoa_ESDM_raw*FC_resample 
plot(FC_habitatSui)

#.................................................................................
## Reclassify raster based on percentiles to delineate habitat suitability classes 
#.................................................................................
# Duplicate distribution model output for editing
zAnoa <- FC_habitatSui

# Define percentile intervals of interest, corresponds to the desired number of habitat suitability classes
# Here, we select percentiles at intervals of 20, which will provide five habitat suitability classes
# Note, percentiles must always start at zero, and omit the max value (i.e. 1)
pp <- global(zAnoa, fun=quantile, na.rm=T,probs=seq(0, 0.8, 0.2))

## Reformat raster for processing (convert from SpatRaster to raster object)
zAnoa_old <- raster::raster(zAnoa)

# Retrieve all raster cell values from distribution model output
gv <- raster::getValues(zAnoa_old)

# Define percentile thresholds based on distribution of values within the SDM layer
ix_A <- findInterval(gv, pp)

# Reclassify distribution model based on threshold values
cAnoa <- setValues(zAnoa_old, ix_A)
plot(cAnoa)

## Reformat back to SpatRaster
cAnoa2 <- rast(cAnoa)
plot(cAnoa2)

# Export the raster layer for postprocessing and land cover calculations in qGIS
writeRaster(cAnoa2, "outputs/anoa_HabSuit.tiff", overwrite=T)

