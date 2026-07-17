# =========================================================
# 3_Complementarity_dimension.R
# =========================================================
#
# DESCRIPTION:
# This script calculates the "complementarity dimension"
# of the Geographic Survey Index (GSI). This dimension
# evaluates the completeness of biological sampling within
# each spatial unit of the study area by comparing observed
# species richness against estimated species richness.
#
# The approach assumes that well-sampled regions should
# contain a number of observed species close to the
# estimated total richness expected for that locality.
# Conversely, regions where estimated richness greatly
# exceeds observed richness are interpreted as
# under-sampled areas.
#
# Species richness is estimated using non-parametric
# richness estimators, specifically Bootstrap and
# Jackknife methods, following the logic that inventory
# completeness can be quantified as:
#
#     Observed richness / Estimated richness
#
# Values closer to 1 indicate more complete sampling,
# whereas lower values indicate sampling gaps.
#
# INPUTS:
# - A shapefile defining the Area of Interest (AOI)
# - A table of biological occurrence records containing:
#     * Record identifiers
#     * Latitude and longitude coordinates
#     * Species names
#
# - Auxiliary functions stored in:
#
#     GAPfunctions.R
#
# OUTPUTS:
# - Raster layers representing sampling completeness
#   estimated using:
#
#     * Bootstrap richness estimator
#     * Jackknife richness estimator
#
# - Output raster files:
#
#     3_Complementarity/Complementariedad_Jacknife.tif
#     3_Complementarity/Complementariedad_Bootstrap.tif
#
# - An RData workspace containing all intermediate objects:
#
#     3_Complementarity/complementariedad_R_object.RData
#
# MAIN PROCESSING STEPS:
#
# 1. Data loading and preprocessing
#    - Load the study area shapefile
#    - Load biological occurrence records
#    - Remove records with missing coordinates
#    - Remove duplicated localities
#
# 2. Spatial data preparation
#    - Convert occurrence records into spatial objects
#    - Rasterize the Area of Interest (AOI)
#    - Assign each occurrence record to a raster cell
#
# 3. Species-by-cell organization
#    - Build species lists for each raster cell
#    - Calculate the number of records per spatial unit
#    - Filter cells according to a minimum record threshold
#
# 4. Richness estimation
#    - Estimate expected species richness using:
#         * Bootstrap estimator
#         * Jackknife estimator
#
# 5. Sampling completeness calculation
#    - Calculate completeness as:
#
#         observed richness / estimated richness
#
#    - Generate completeness values for each raster cell
#    - Constrain values greater than 1 to a maximum of 1
#
# 6. Raster generation and visualization
#    - Convert completeness estimates into raster layers
#    - Generate histograms and density distributions of
#      completeness values
#
# 7. Export results
#    - Save raster outputs for Bootstrap and Jackknife
#      completeness layers
#    - Save the R workspace for reproducibility
#
# =========================================================

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
Data <- Data[, c('gbifID', 'decimalLatitude', 'decimalLongitude', 'species')] # Seleccionar las comlumnas donde la 1ra sea el id del registro, 2da y 3ra sea la latitud y longitud y la 4ta el nombre de las especies.
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
spListByCell <- Data[!is.na(Data$celdas), c('species', 'celdas')]
spListByCell <- na.omit(spListByCell) # Eliminar NAs
freqTable <- table(spListByCell$celdas) # Contar frecuencia de especies por celda

# Definir umbral de frecuencia
treshold <- 0

# Filtrar celdas que cumplen con el umbral de frecuencia
spListByCell <- spListByCell[spListByCell$celdas %in% names(which(freqTable >= treshold)), ]

# Estimar la riqueza usando los métodos de bootstrap y jackknife
estimateS <- richEst(sppList = spListByCell$species, indexID = spListByCell$celdas)

rm(spListByCell)

# Crear rasters para las estimaciones
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
