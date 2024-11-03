#rm(list = ls(all = TRUE))

library(terra)
library(sf)
library(janitor)
library(dplyr)
library(vroom)

# ---------------------------
#   Data loading
# ---------------------------
# Establecer directorio de trabajo y crear directorio para resultados
wd <- ''
setwd(wd)                              # Cambiar el directorio de trabajo a la ruta especificada
dir.create("3_Complementarity/")       # Crear un nuevo directorio para guardar los resultados

# Cargar funciones adicionales desde un archivo externo
source('./GAPfunctions.R')   # Cargar funciones desde el archivo 'GAPfunctions.R'

# Definir el sistema de referencia de coordenadas (CRS)
GRS.geo <- '+proj=longlat +datum=WGS84 +no_defs'

#col<-getData(name = 'GADM', country = 'COL', level = 0) #Areo of interes (AOI) # función desactivada. Sirve para cargar shapefile de Colombia desde servidor remoto del paquete geodata
col <- read_sf('/Carpeta_shapefile', 'shapefile')
colrast <- rast(col, res=0.008333333) # Crear un raster con la resolución especificada
colrast <- terra::rasterize(col, colrast, fun="sum") # Rasterizar el shapefile

# Carga de registros
Data <- vroom(" ", col_names = TRUE)
Data <- Data[, c('gbifID', 'decimalLatitude', 'decimalLongitude', 'species')] # Seleccionar las comlumnas donde la 1ra sea el id del registro, 2da t 3ra sea la latitud y longitud y la 4ta el nombre de las especie.
Data <- Data[!is.na(Data$decimalLongitude_x),] # Eliminar registros con longitud faltante
Data <- Data[!is.na(Data$decimalLatitude_x),]  # Eliminar registros con latitud faltante
Data2 = Data %>% distinct(decimalLatitude_x, decimalLongitude_x, .keep_all = TRUE) # Eliminar duplicados por coordenadas

# Convertir a objeto 'sf'
Data2 <- st_as_sf(Data2, coords = c('decimalLongitude', 'decimalLatitude'), crs = GRS.geo)

# Crear una grilla con el raster y extraer celdas para los puntos
grilla <- colrast
grilla[] <- 1:ncell(colrast) # Asignar un identificador único a cada celda del raster
en_area <- mask(grilla, vect(col)) # Aplicar una máscara para limitar al área de interés
celdas <- terra::extract(en_area, vect(Data2)) # Extraer identificadores de celda para cada punto
celdas = celdas[,"layer"] # Seleccionar la columna de celdas
Data2$celdas <- celdas # Añadir columna de celdas a los datos
Data <- as.data.frame(Data2) # Convertir a dataframe
rm(celdas) # Eliminar objeto 'celdas' para liberar memoria

# ---------------------------
# 2. Data base completness
# ---------------------------

# Crear una tabla de frecuencia para las celdas
spListByCell <- Data[!is.na(Data$celdas), c('species_x', 'celdas')]
spListByCell <- na.omit(spListByCell) # Eliminar NA
freqTable <- table(spListByCell$celdas) # Contar frecuencia de especies por celda

# Definir umbral de frecuencia
treshold <- 0

# Filtrar celdas que cumplen con el umbral de frecuencia
spListByCell <- spListByCell[spListByCell$celdas %in% names(which(freqTable >= treshold)), ]

# Estimar la riqueza usando los métodos de bootstrap y jackknife
estimateS <- richEst(sppList = spListByCell$species, indexID = spListByCell$celdas)

rm(spListByCell)

# Inicializar rasters para las estimaciones
compRichBoot <- compRichJack <- richJackHQ <- richBootHQ <- richJack <- richBoot <- en_area * 0
richBoot[as.numeric(rownames(estimateS))] <- estimateS$Boot # Asignar valores de bootstrap
richJack[as.numeric(rownames(estimateS))] <- estimateS$JNhat # Asignar valores de jackknife

compRichBoot[as.numeric(rownames(estimateS))] <- estimateS$Sobs/estimateS$Boot # Calcular la razón de riqueza bootstrap
compRichJack[as.numeric(rownames(estimateS))] <- estimateS$Sobs/estimateS$JNhat # Calcular la razón de riqueza jackknife

# Ajustar valores de complejidad a 1 o NA según la presencia de datos
#compRichBoot[as.numeric(rownames(estimateS)[estimateS$Sobs == 1])] <- NA
#compRichJack[as.numeric(rownames(estimateS)[estimateS$Sobs == 1])] <- NA
compRichJack[compRichJack[] >= 1] <- 1
compRichBoot[compRichBoot[] >= 1] <- 1

# Filtrar valores válidos y generar histogramas
richBootVals <- compRichBoot[!is.na(compRichBoot[]) & compRichBoot[] != 0]
richJackVals <- compRichJack[!is.na(compRichJack[]) & compRichJack[] != 0]
hist(richBootVals, main = 'Density Bootstrap', freq = FALSE, xlim = c(0, 1.2))
lines(density(richBootVals), main = 'Density Bootstrap')

# Guardar resultados
writeRaster(compRichBoot, paste0("3_Complementarity/Complementariedad_Jacknife.tif"), overwrite=TRUE)
writeRaster(compRichJack, paste0("3_Complementarity/Complementariedad_Bootstrap.tif"), overwrite=TRUE)

# Guardar el entorno de trabajo
save.image(paste0('3_Complementarity/complementariedad_R_object.RData'))
