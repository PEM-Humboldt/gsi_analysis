
library(sf)
library(maptools)
library(spatstat)
library(raster)

# ---------------------------
#   Data loading
# ---------------------------
wd <- 'F:/GSI/'
setwd(wd)
dir.create("1_Records/")

col<-getData(name = 'GADM', country = 'COL', level = 0) #Areo of interes (AOI)
col@data$OBJECTID<-1
GRS.geo<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Geographic Reference to AOI
CRS.proj<-'+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' # Planar Reference to AOI
proj.col<- spTransform(col, CRS.proj)
shape_zoneOwin<-as.owin(proj.col)

load('data_to_gap.RData') # Records to use in the analysis / require name of the species, coordinates (geograhpic type) and ID identifier
Data <- data_to_gap
Data2 <- Data[, c('ID', 'lat', 'lon')]
Data2 <- Data2[!is.na(Data2$lat),]
Data2 <- Data2[!is.na(Data2$lon),]
coordinates(Data2)=~)on+lat

coords<-as.data.frame(unique(Data2@coords))
coords$ID<-1:nrow(coords)
coordinates(coords)=~lon+lat
coords@proj4string@projargs<-GRS.geo
coordinates.col<- spTransform(coords, CRS.proj)
system.time(over.coords <- over(coordinates.col, proj.col)); head(over.coords) #Identify coordinates within the AOI

enAOI <- as.numeric(over.coords$OBJECTID)
enAOI[enAOI>=1] <- 1
enAOI[is.na(enAOI)] <- 0
coordinates.col<-as.data.frame(coordinates.col)
coordinates.col$enAOI <- enAOI
coordinates.col<- coordinates.col[coordinates.col$enAOI == 1,]


#shape_zoneOwin <- as.owin(col)
p <- ppp(coordinates.col$lon, coordinates.col$lat, window=shape_zoneOwin, unitname=c("metre","metres")) # Create a Point Pattern using unit to GRS
summary(p)
plot(p)

diggle <- bw.diggle(p) # Selected value 
plot(diggle)
plot(diggle, xlim= c(0,100), main="Smoothing bandwidth for the kernel estimation")

system.time(diggle_den <- density.ppp(p, diggle, eps=1000)) # Kernel Smoothed Intensity of Point Pattern 1km - 17 min
plot(diggle_den,diggle, main='Gap density diggle 1km')

densi_dig1km<-raster(diggle_den)
rescal_dig1km<-densi_dig1km/max(diggle_den) # Re-scaling the value using teh Max value. Please, change the value based on your results.

crs(rescal_dig1km)<-CRS.proj
rescal_dig1km_wgs84<- projectRaster(rescal_dig1km, crs = GRS.geo)
plot(rescal_dig1km_wgs84)
writeRaster(rescal_dig1km_wgs84, filename="1_Records/Vac_dens_rescal_1km.tif", format="GTiff", overwrite=TRUE)
