# =========================================================
# 1_Record_dimension.R
# =========================================================
#
# DESCRIPTION:
# This script calculates the "record dimension" of the
# Geographic Survey Index (GSI), following the framework
# proposed by García Márquez et al. (2012). This dimension
# represents the spatial concentration of biological records
# and is used as a proxy for sampling effort across the
# study region.
#
# The workflow uses georeferenced occurrence records to
# estimate a kernel density surface, where areas with a
# high concentration of records are interpreted as highly
# sampled regions, while areas with low record density are
# interpreted as poorly sampled regions.
#
# INPUTS:
# - A shapefile defining the Area of Interest (AOI)
# - A table of biological records containing:
#     * Unique record identifiers
#     * Latitude coordinates
#     * Longitude coordinates
#
# OUTPUTS:
# - A raster layer representing the normalized density of
#   biological records across the study area:
#
#     1_Records/Vac_dens_rescal_1km.tif
#
# - An RData file containing the workspace generated during
#   the analysis:
#
#     1_Records/Records_R_object.RData
#
# MAIN PROCESSING STEPS:
#
# 1. Data loading and preprocessing
#    - Load the study area shapefile
#    - Load biological occurrence records
#    - Remove records with missing coordinates
#    - Remove duplicated geographic coordinates
#
# 2. Coordinate system standardization
#    - Transform all spatial data into a projected
#      coordinate system suitable for spatial analysis
#
# 3. Spatial filtering
#    - Retain only records located within the Area of
#      Interest (AOI)
#
# 4. Point pattern construction
#    - Convert occurrence coordinates into a spatial
#      point pattern object compatible with spatstat
#
# 5. Kernel density estimation
#    - Estimate sampling intensity using kernel density
#      smoothing
#    - Select the smoothing bandwidth using the Diggle
#      method
#
# 6. Raster generation and normalization
#    - Convert the kernel density output into a raster
#    - Rescale raster values between 0 and 1
#    - Reproject raster into WGS84 geographic coordinates
#
# 7. Export results
#    - Save the final raster layer
#    - Save the R workspace for reproducibility
#
# =========================================================


library(sf)
library(terra)
library(spatstat)
library(geodata)
library(dplyr)
library(vroom)


# ---------------------------
#   Data loading
# ---------------------------

# Definición del directorio de trabajo
wd <- '/GSI'
setwd(wd)  # Establecer el directorio de trabajo
dir.create("1_Records/")  # Crear un directorio para guardar resultados

# Cargar shapefile del área de estudio
##col <- geodata::gadm(country= "COL", level=0, path= ".") # función desactivada. Sirve para cargar shapefile de Colombia desde servidor remoto del paquete geodata
col <- read_sf('/Carpeta_shapefile', 'shapefile')

# Definición del sistema de referencia geográfico (GCS) y proyectado (CRS)
GRS.geo<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Geographic Reference to AOI
CRS.proj<-'+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' # Planar Reference to AOI

# Transformación del shapefile al sistema de coordenadas proyectadas
proj.col <- sf::st_transform(col, crs = CRS.proj) 

# Convertir el shapefile a una ventana espacial compatible con spatstat
shape_zoneOwin <- as.owin(proj.col)

# Cargar registros de especies desde un archivo txt
Data <- vroom("Archivo_registros.txt", col_names = TRUE)
Data2 <- Data[, c('gbifID', 'decimalLatitude', 'decimalLongitude')] #Select ad first column the name of the id record. And for second and third column, latitude and longitude, respectively.
Data2 <- Data2[!is.na(Data2$decimalLatitude),] # Eliminar registros con latitud NA
Data2 <- Data2[!is.na(Data2$decimalLongitude),] # Eliminar registros con longitud NA

# Eliminar duplicados basados en coordenadas
Data2 <- Data2[!duplicated(Data2[c("decimalLongitude", "decimalLatitude")]), ]
# De acuerdo con García Márquez et al. (2012), esta dimensión se basa en la lista de localidades de colecta como puntos
# con los cuales se genera la "densidad de localidades" ("density of collection localities").

# Convertir a objeto sf y transformar las coordenadas al sistema proyectado
coords <- st_as_sf(Data2, coords = c("decimalLongitude", "decimalLatitude"), crs = GRS.geo)
coordinates.col <- st_transform(coords, crs = CRS.proj)

# Intersección espacial para verificar si las coordenadas están dentro del AOI (Area of Interest)
system.time(over.coords <- st_intersects(coordinates.col, proj.col, sparse = FALSE)); head(over.coords)
Data2 <- Data2[over.coords == TRUE,]
#system.time(over.coords <- st_intersection(coordinates.col, proj.col, sparse = FALSE)); head(over.coords) #Identify coordinates within the AOI
#coordinates.col <- over.coords #Este paso asume que el st_intersect solo genera en su output coordenadas que efextivamente haver overlap (no genera NAs).
# El comando con st_intersection genera el objetvo espacial de los overlap, pero toma más tiempo.
#system.time(over.coords <- st_within(coordinates.col, proj.col)) #otra alternativa, mirar cuál es más rápida


# Crear un patrón de puntos con la ventana espacial definida
coords_df <- st_coordinates(coordinates.col)
p <- ppp(coords_df[, "X"], coords_df[, "Y"], window = shape_zoneOwin, unitname = c("metre", "metres"))
summary(p)
plot(p)

# Estimación de la densidad del kernel
diggle <- bw.diggle(p)  # Selección del valor de suavizamiento
plot(diggle)
plot(diggle, xlim= c(0,100), main="Smoothing bandwidth for the kernel estimation")

# Calcular la intensidad suavizada del patrón de puntos (1km de resolución)
system.time(diggle_den <- density.ppp(p, diggle, eps=10000)) # Kernel Smoothed Intensity of Point Pattern 1km - 17 min

# Graficar la densidad suavizada
plot(diggle_den,diggle, main='Gap density diggle 1km')

# Convertir la densidad a un objeto raster
densi_dig1km <- rast(diggle_den)
rescal_dig1km <- densi_dig1km/max(diggle_den)  # Reescalar usando el valor máximo

# Asignar el sistema de referencia proyectado al raster
crs(rescal_dig1km) <- CRS.proj

# Proyectar el raster reescalado al sistema de referencia geográfico (WGS84)
rescal_dig1km_wgs84 <- terra::project(rescal_dig1km, crs(GRS.geo))

# Graficar el raster proyectado
plot(rescal_dig1km_wgs84)

# Guardar el raster resultante
writeRaster(rescal_dig1km_wgs84, filename="1_Records/Vac_dens_rescal_1km.tif", overwrite=TRUE)

# Guardar el entorno de trabajo en un archivo RData
save.image(paste0('1_Records/Records_R_object.RData'))
