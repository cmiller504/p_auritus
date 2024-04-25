## Generalized Dissimilarity modelling - data preparation, model fitting, testing  
## P.auritus SNP data with FST table                                                

library(gdm)

# Load the environmental grids
rasters <- c("EVs/bio01.asc",
             "EVs/bio04.asc",
             "EVs/bio12.asc",
             "EVs/bio15.asc",
             "EVs/bio19.asc")

# stack variables
raster_stack <- raster::stack(rasters)

# set up resistance matrices
resistDist = list()
resistDist[[1]] = read.csv("inputs/resistance_matrices/ElevationRivers100.csv")
resistDist[[2]] = read.csv("inputs/resistance_matrices/ElevationRiversScaled.csv")
resistDist[[3]] = read.csv("inputs/resistance_matrices/miroc.csv")
resistDist[[4]] = read.csv("inputs/resistance_matrices/ccsm.csv")


# Read in the input table with FST values
gdmTab_FST <- read.csv("GDM_input_table_maf_rSNP_fst.csv")


# Fit a GDM on the input table we just created with geographic distances included

gdm.1.0<-gdm(gdmTab_FST, geo=TRUE)

#summary of the GDM
summary(gdm.1.0)


# Plot model showing functional response curves for significant environmental predictors
plot(gdm.1.0, 
     include.rug=TRUE, 
     rug.sitepair=gdmTab_FST)

## Test the significance of the model & its variables *************************

gdm.1.0.test<- gdm.varImp(gdmTab_FST, geo=TRUE, nPerm = 1000, parallel = TRUE)


## Validate the GDM 

fst_crossval =gdm.crossvalidation(gdmTab_FST, train.proportion=0.7, n.crossvalid.tests=100,geo=TRUE, splines=NULL, knots=NULL)
