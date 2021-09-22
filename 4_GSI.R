
library(raster)
library(rgdal)
library(hyperSpec)
# ---------------------------
#   Data loading
# ---------------------------
wd <- 'F:/GSI/'
setwd(wd)
dir.create("GSI/")


GRS.geo<- '+proj=longlat +datum=WGS84 +no_defs'
col<-getData(name = 'GADM', country = 'COL', level = 0) #Areo of interes (AOI)
col@data$OBJECTID<-1
colr<-raster(extent(col), res=0.008333333)
colrast<- rasterize(col,colr)

# ---------------------------
# 4. GSI / Gap selection index
# ---------------------------

# Estandarizacion de valores

#densRegistros <- raster('INDICE_DENSIDAD.tif')
#d <- raster('INDICE_AMBIENTAL.tif')
DENSIDAD <- raster("1_Records/Vac_dens_rescal_1km.tif")
crs(DENSIDAD)<-GRS.geo
AMBIENTAL <- raster("2_Ambiental/ambiental_dimension_GSI_rescal.tif")
crs(AMBIENTAL)<-GRS.geo

COMPLEMENTARIEDAD_BOOT<-raster("3_Complementarity/Complementariedad_Bootstrap.tif")
crs(COMPLEMENTARIEDAD_BOOT)<-GRS.geo
COMPLEMENTARIEDAD_JACK<-raster("3_Complementarity/Complementariedad_Jacknife.tif")
crs(COMPLEMENTARIEDAD_JACK)<-GRS.geo


ambiental_resample <-resample(AMBIENTAL, DENSIDAD)
comp_jack <-resample(COMPLEMENTARIEDAD_JACK, DENSIDAD)
comp_boot <-resample(COMPLEMENTARIEDAD_BOOT, DENSIDAD)

AMBValsNoNA <- AMBIENTAL[!is.na(AMBIENTAL[])]
AMBValsNoNANoZeros <- AMBIENTAL[!is.na(AMBIENTAL[]) & AMBIENTAL[] != 0]
hist(AMBValsNoNA, main = 'Ambiental no NA')
hist(AMBValsNoNANoZeros, main = 'Ambiental /n no NA no 0')

# # Completness
par(mfrow = c(1, 2))
plot(COMPLEMENTARIEDAD_BOOT, main = 'Bootstrap')
plot(col, add = TRUE)
plot(COMPLEMENTARIEDAD_JACK, main = 'Jackknife')
plot(col, add = TRUE)

#  GSI / GAP INDEX
rasterOptions(tolerance = 0.1)
GSI_BOOT <- (3- DENSIDAD - ambiental_resample- comp_boot)/3
GSI_JACK <- (3-DENSIDAD - ambiental_resample - comp_jack)/3

par(mfrow = c(1, 2))

plot(GSI_BOOT, main = 'GAP SELCTION INDEX (BOOT)', zlim = c(0, 1))
plot(col, add = TRUE)

plot(GSI_JACK, main = 'GAP SELCTION INDEX (JACK)', zlim = c(0, 1))
plot(col, add = TRUE)

densJack <- density(GSI_JACK, main = 'Jackk', xlim = c(0, 1))
densBoot <- density(GSI_BOOT, main = 'Boot', xlim = c(0, 1))

par(mfrow = c(1, 1))
plot(densJack, main = 'GSI values', col = 'blue', ylim = c(0, 10))
lines(densBoot, col = 'red')
legend('topleft', legend = c('Jackknife', 'Bootstrap'), 
       lty = c(1, 1), lwd = c(1, 1), col = c('blue', 2))


writeRaster(GSI_BOOT,"GSI/INDICE_GSI_BOOT.tif", overwrite = TRUE)
writeRaster(GSI_JACK,"GSI/INDICE_GSI_JACK.tif", overwrite = TRUE)
