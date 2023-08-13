## R Studio code to support Buckland et al. 2023 GCB Bioenergy manuscript

## Please see below code for example of how each species was projected using the biomod2 package - updates according to package changes; network mapping; and different species will need to be changed by the individual user
## The code below is largely adapted from the examples given on the biomod2 package website and documentation

## biomod2 citation: Thuiller W, Georges D, Gueguen M, Engler R, Breiner F, Lafourcade B, Patin R. 2023. biomod2: Ensemble Platform for Species Distribution Modeling. R package version 4.2-4-5.

# 1. Load required libraries
library(biomod2)
library(dismo)
library(rJava)
library(raster)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)
library(countrycode)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(rlang)

# 2. Obtain occurrence data from GBIF via rgbif
z <- occ_data(scientificName = "Opuntia ficus-indica", basisOfRecord = 'HUMAN_OBSERVATION', limit = 10000,
              hasCoordinate = T)
dat <- data.frame(z$data)

# 3. Select columns of interest from occurrence record
dat <- dat %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)

# 4. Remove records without coordinates
dat <- dat%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# 5. Plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray75")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.5)+
  theme_bw()

# 6. Convert country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

# 7. Flag problems with occurrence data using clean_coordinates package
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
                           countries = "countryCode", 
                           species = "species",
                           tests = c("capitals", "duplicates", "seas", "centroids", "equal","gbif", "institutions",
                                     "zeros"))
summary(flags)
plot(flags, lon = "decimalLongitude", lat="decimalLatitude")

dat_cl <- dat[flags$.summary,]

wm <- borders("world", colour="gray50", fill="gray75")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = flags, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "black", size = 0.5)+
  geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "yellow", size = 0.5)+
  theme_bw()

head(dat_cl)

data_subset <- dat_cl[,c(2,3)]

data_subset$OpFI <-rep(1,nrow(data_subset))

head(data_subset)

# 8. Name occurrence data columns 
names(data_subset)[1] <- "X_WGS84"
names(data_subset)[2] <- "Y_WGS84"
head(data_subset)

myRespName <- 'OpFI_test'

# 9. Save presence data for our species to a variable
myResp <- as.numeric(data_subset[,c(3)])

# 10. XY coordinates of species data
myRespXY <- data_subset[,c("X_WGS84","Y_WGS84")]

# 11. Set extent of occurrence data
OFIext = extent(-180, 180, -90, 90)
r <- raster(OFIext)
res(r) <- 1

# 12. Sub-sample occurrence data to a maximum of points per grid cell. Example n=5
myRespXYsel <- gridSample(myRespXY,r,n=5)

rownum <- nrow(myRespXYsel)

plot(myRespXY)
points(myRespXYsel,cex=1,col='red',pch='x')

wm <- borders("world", colour="gray50", fill="gray75")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = myRespXYsel, aes(x = X_WGS84, y = Y_WGS84),
             colour = "darkgreen", size = 0.5)+
  theme_bw()

myRespCut <-head(myResp,n=rownum)

# 13. Load required bioclim rasters and stack as a predictor variable dataset
WorldClimData <- getData('worldclim',var='bio',res=2.5)
bio_02 <- subset(WorldClimData,2)
bio_06 <- subset(WorldClimData,6)
bio_12 <- subset(WorldClimData,12)
bio_15 <- subset(WorldClimData,15)
predictors <- stack(bio_02, bio_06, bio_12, bio_15)

# 14. Crop to Africa extent for projecting and create cropped extent predictor raster stack
ext = extent(-22, 52, -38, 38)
bio_02c <- crop(bio_02,ext)
bio_06c <- crop(bio_06,ext)
bio_12c <- crop(bio_12,ext)
bio_15c <- crop(bio_15,ext)
predictors_cut <- stack(bio_02c, bio_06c, bio_12c, bio_15c)

# 15. Build SDM based on historical/current climate

# 16. Format the data and add background pseudo-absence data points
myBiomodData <- BIOMOD_FormatingData(resp.var = myRespCut,
                                     expl.var = predictors,
                                     resp.xy = myRespXYsel,
                                     resp.name = myRespName,
                                     eval.resp.var = NULL,
                                     eval.expl.var = NULL,
                                     eval.resp.xy = NULL,
                                     PA.nb.rep = 5,
                                     PA.nb.absences = rownum,
                                     PA.strategy='sre',
                                     PA.dist.min = 20000,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     #PA.table = NULL,
                                     na.rm = TRUE)
plot(myBiomodData)

# 17. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips=list(path_to_maxent.jar = 'INSERT PATH TO MAXENT FILE HERE'))  

# 18. Computing the models
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = c('GBM', 'RF'),
  bm.options = myBiomodOption,
  nb.rep = 10,
  data.split.perc =60,
  prevalence=NULL,
  var.import=3,
  metric.eval = c('TSS','ROC'),
  save.output = TRUE,
  #compress = 'xz',
  #rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

myBiomodModelOut

# 19. Get all models evaluation and print dimension names of object
myBiomodModelEval_current <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval_current)

# 20. Calculate average variable importance metrics by algorithm and then by variable
variables_imp <- get_variables_importance(myBiomodModelOut)
dimnames(variables_imp)

variables_imp[,"GBM",,]
mean(variables_imp["bio_02","GBM",,])
mean(variables_imp["bio_06","GBM",,])
mean(variables_imp["bio_12","GBM",,])
mean(variables_imp["bio_15","GBM",,])

variables_imp[,"RF",,]
mean(variables_imp["bio_02","RF",,])
mean(variables_imp["bio_06","RF",,])
mean(variables_imp["bio_12","RF",,])
mean(variables_imp["bio_15","RF",,])

mean(variables_imp["bio2",,])
mean(variables_imp["bio6",,])
mean(variables_imp["bio12",,])
mean(variables_imp["bio15",,])

# 21. Build ensemble model and view evaluation metrics
myBiomodEM_current_all <- BIOMOD_EnsembleModeling(
  myBiomodModelOut,
  models.chosen = 'all',
  em.by='all',
  metric.select = 'all',
  metric.select.thresh = c(0.4,0.7),
  metric.eval = c('TSS','ROC'),
  prob.mean = F,
  prob.mean.weight = T,
  prob.cv = T,
  prob.ci = F,
  prob.ci.alpha = 0.05,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_current_all

get_evaluations(myBiomodEM_current_all)
get_built_models(myBiomodEM_current_all)
get_kept_models(myBiomodEM_current_all, model=1)

# 22. Projection over the cropped extent under current conditions
myBiomodProj_current <- BIOMOD_Projection(
  myBiomodModelOut,
  new.env = predictors_cut,
  proj.name = 'Current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

# 23. Make plots sub-selected by str.grep argument
plot(myBiomodProj_current, str.grep = 'GBM')
plot(myBiomodProj_current, str.grep = 'RF')

# 24. Creating the ensemble projections
myBiomodEF_current <- BIOMOD_EnsembleForecasting( 
  EM.output = myBiomodEM_current_all, 
  projection.output = myBiomodProj_current,
  binary.meth = 'TSS',
  output.format = ".img")

plot(myBiomodEF_current)

# 25. Create example 2D and 3D response curves 
rp.dat.curves <- bm_PlotResponseCurves(
  bm.out = myBiomodEM_current_all,
  models.chosen = 'OpFI.test_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData',
  show.variables = c('bio2','bio6','bio12','bio15'),
  do.bivariate = FALSE,
  fixed.var = 'mean')

rp.dat.curves.3d <- bm_PlotResponseCurves(
  bm.out = myBiomodEM_current_all,
  models.chosen = 'OpFI.2023b_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData',
  show.variables = c('bio2','bio6','bio12','bio15'),
  do.bivariate = TRUE,
  fixed.var = 'mean')

## Second section of analysis then takes the models trained based on current predictor-occurrence relationships and re-projects SDM outputs based on new environmental predictor data stacks 

# Please see below example for one SSP scenario - user will need to insert details to where relevant future bioclim projections are stored on their device

# 26. Load new dataset SSP126 2081-2100
bio_02 <- raster("INSERT FILE LOCATION")
bio_06 <- raster("INSERT FILE LOCATION")
bio_12 <- raster("INSERT FILE LOCATION")
bio_15 <- raster("INSERT FILE LOCATION")

# 27. Create FULL EXTENT predictor stacks
predictors_SSP126_2081_2100 <- stack(bio_02, bio_06, bio_12, bio_15)

# 28. Crop to Africa extent for projecting
ext = extent(-22, 52, -38, 38)
bio_02c <- crop(bio_02,ext)
bio_06c <- crop(bio_06,ext)
bio_12c <- crop(bio_12,ext)
bio_15c <- crop(bio_15,ext)

predictors_SSP126_2081_2100_cut <- stack(bio_02c, bio_06c, bio_12c, bio_15c)

# 29. Projection over the globe under new scenario conditions using trained models from earlier steps
myBiomodProj_SSP126_2081_2100 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = predictors_SSP126_2081_2100_cut,
  proj.name = 'SPP126_2081_2100',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

# 30. Creating the future SSP ensemble projections
myBiomodEF_SSP126_2081_2100 <- BIOMOD_EnsembleForecasting( 
  EM.output = myBiomodEM_current, 
  projection.output = myBiomodProj_SSP126_2081_2100,
  binary.meth = 'TSS',
  output.format = ".img")

plot(myBiomodEF_SSP126_2081_2100)