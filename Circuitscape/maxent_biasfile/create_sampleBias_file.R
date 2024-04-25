### Create a sampling bias file to use in Maxent SDM

library(dismo) # interface with MaxEnt
library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools) # reading shapefiles

#load occurrence data

occurdat = read.csv("Pauritus_GPS.csv")

locations = cbind(occurdat[,3:2])

# need occurrences on a spatial grid
# turn them into a raster

climdat <- raster("elev.asc")
crs(climdat) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

occur.ras <- rasterize(locations, climdat, field = 1)
plot(occur.ras)

# mask to shapefile of land

e <- as(extent(8, 17, -5, 14), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
land = shapefile("ne_10m_land/ne_10m_land.shp")
mask = crop(land,e)
  
occur.mask <- mask(occur.ras, mask) %>% crop(mask)

# make two-dimensional kernel density estimate to get bias file

presences <- which(values(occur.mask) == 1)
pres.locs <- coordinates(occur.mask)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.mask), ncol(occur.mask)))
dens.ras <- raster(dens)
plot(dens.ras)

writeRaster(dens.ras, "sample_bias.tif")
writeRaster(dens.ras, "sample_bias.grd")
