#rm(list = ls(all = TRUE))
library(raster)
library(rgdal)
library(janitor)
library(dplyr)

# ---------------------------
#   Data loading
# ---------------------------
wd <- 'F:/GSI/'
setwd(wd)
dir.create("3_Complementarity/")

source('GAPfunctions.R')

GRS.geo<- '+proj=longlat +datum=WGS84 +no_defs'
col<-getData(name = 'GADM', country = 'COL', level = 0) #Areo of interes (AOI)
col@data$OBJECTID<-1
colr<-raster(extent(col), res=0.008333333)
colrast<- rasterize(col,colr)

#Records
load('data_to_gap.RData')
Data <- data_to_gap
Data <- Data[, c('ID', 'lat', 'lon', 'species')]
Data <- Data[!is.na(Data$lon),]
Data <- Data[!is.na(Data$lat),]
Data$coords <- paste0(Data$lat, ' ', Data$lon)
Data2<- Data[!duplicated(Data$coords),]
coordinates(Data2)=~lon+lat
Data2@proj4string@projargs <- GRS.geo

grilla <- colrast
grilla[] <- 1:ncell(colrast)
en_area <- mask(grilla, col)
celdas <- extract(en_area, Data2)
Data2$celdas <- celdas
Data<-as.data.frame(Data2)
rm(celdas)

# ---------------------------
# 2. Data base completness
# ---------------------------

spListByCell <- Data[!is.na(Data$celdas), c('origin', 'celdas')]
freqTable <- table(spListByCell$celdas)

treshold <- 0

spListByCell <- spListByCell[spListByCell$celdas %in% names(which(freqTable >= treshold)), ]
estimateS <- richEst(sppList = spListByCell$species, indexID = spListByCell$celdas)

rm(spListByCell)

compRichBoot <- compRichJack <- richJackHQ <- richBootHQ <- richJack <- richBoot <- en_area * 0
richBoot[as.numeric(rownames(estimateS))] <- estimateS$Boot
richJack[as.numeric(rownames(estimateS))] <- estimateS$JNhat

compRichBoot[as.numeric(rownames(estimateS))] <- estimateS$Sobs/estimateS$Boot
compRichJack[as.numeric(rownames(estimateS))] <- estimateS$Sobs/estimateS$JNhat

#compRichBoot[as.numeric(rownames(estimateS)[estimateS$Sobs == 1])] <- NA
#compRichJack[as.numeric(rownames(estimateS)[estimateS$Sobs == 1])] <- NA
compRichJack[compRichJack[] >= 1] <- 1
compRichBoot[compRichBoot[] >= 1] <- 1

richBootVals <- compRichBoot[!is.na(compRichBoot[]) & compRichBoot[] != 0]
richJackVals <- compRichJack[!is.na(compRichJack[]) & compRichJack[] != 0]
hist(richBootVals, main = 'Density Bootstrap', freq = FALSE, xlim = c(0, 1.2))
lines(density(richBootVals), main = 'Density Bootstrap')


writeRaster(compRichBoot, paste0("3_Complementarity/Complementariedad_Jacknife.tif"), overwrite=TRUE)
writeRaster(compRichJack, paste0("3_Complementarity/Complementariedad_Bootstrap.tif"), overwrite=TRUE)

save.image(paste0('3_Complementarity/complementariedad_R_object.RData'))
