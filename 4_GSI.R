
library(terra)
library(sf)
library(hyperSpec)
library(RColorBrewer)

# ---------------------------
#   Data loading
# ---------------------------
# Establecer el directorio de trabajo y crear un directorio para resultados
wd <- ''
setwd(wd)                              # Cambiar el directorio de trabajo a la ruta especificada
dir.create("GSI/")                     # Crear un nuevo directorio para guardar los resultados

# Definir el sistema de referencia de coordenadas (CRS)
GRS.geo <- '+proj=longlat +datum=WGS84 +no_defs'

# Leer el shapefile del área de interés (AOI)
###col <- geodata::gadm(country= "COL", level=0, path= ".") #Areo of interes (AOI) # función desactivada. Sirve para cargar shapefile de Colombia desde servidor remoto del paquete geodata
col <- read_sf('/Carpeta_shapefile', 'shapefile')
col <- sf::st_transform(col, GRS.geo)  # Transformar el CRS del shapefile
colrast <- rast(col, res=0.008333333)  # Crear un raster con la resolución especificada
colrast <- terra::rasterize(col, colrast, fun="sum") # Rasterizar el shapefile

# ---------------------------
# 4. GSI / Gap selection index
# ---------------------------

# Estandarización de valores


# Leer los rasters de datos
DENSIDAD <- rast("1_Records/Vac_dens_rescal_1km.tif") # Raster de densidad de registros
crs(DENSIDAD) <- GRS.geo                      # Asignar el CRS al raster de densidad
AMBIENTAL <- rast("2_Ambiental/ambiental_dimension_GSI_rescal_ajusvar.tif") # Raster ambiental
crs(AMBIENTAL) <- GRS.geo                    # Asignar el CRS al raster ambiental

# Leer los rasters de complementariedad
COMPLEMENTARIEDAD_BOOT <- rast("3_Complementarity/Complementariedad_Bootstrap.tif")
crs(COMPLEMENTARIEDAD_BOOT) <- GRS.geo       # Asignar el CRS al raster de bootstrap
COMPLEMENTARIEDAD_JACK <- rast("3_Complementarity/Complementariedad_Jacknife.tif")
crs(COMPLEMENTARIEDAD_JACK) <- GRS.geo        # Asignar el CRS al raster de jackknife

# Reescalar rasters a la misma resolución que el raster de densidad
ambiental_resample <-resample(AMBIENTAL, DENSIDAD)
comp_jack <-resample(COMPLEMENTARIEDAD_JACK, DENSIDAD)
comp_boot <-resample(COMPLEMENTARIEDAD_BOOT, DENSIDAD)

# Filtrar valores no NA y no cero del raster ambiental
AMBValsNoNA <- AMBIENTAL[!is.na(AMBIENTAL[])]
AMBValsNoNANoZeros <- AMBIENTAL[!is.na(AMBIENTAL[]) & AMBIENTAL[] != 0]
hist(AMBValsNoNA, main = 'Ambiental no NA')
hist(AMBValsNoNANoZeros, main = 'Ambiental /n no NA no 0')

# Convertir el objeto 'col' a tipo 'vect' para su uso en gráficos
col=vect(col)

# Graficar las capas de complementariedad
par(mfrow = c(1, 2))
plot(COMPLEMENTARIEDAD_BOOT, main = 'Bootstrap')
lines(col, lwd=1)
plot(COMPLEMENTARIEDAD_JACK, main = 'Jackknife')
lines(col, lwd=1)

# Calcular el índice de selección de huecos (GSI) para ambos métodos
terraOptions(tolerance = 0.1) # Ajustar opciones de tolerancia
GSI_BOOT <- (3 - DENSIDAD - ambiental_resample - comp_boot) / 3 # Índice GSI para bootstrap
GSI_JACK <- (3 - DENSIDAD - ambiental_resample - comp_jack) / 3 # Índice GSI para jackknife


# Definir una paleta de colores y establecer el rango de valores para la visualización
color_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
zlim <- c(0, 1)

# Graficar los índices GSI
par(mfrow = c(1, 2))
plot(GSI_BOOT, main = 'GAP SELCTION INDEX (BOOT)', range = zlim)
lines(col, lwd=1)
plot(GSI_JACK, main = 'GAP SELCTION INDEX (JACK)', range = zlim)
lines(col, lwd=1)

# Graficar densidades de los índices GSI
par(mfrow = c(1, 2))
density(GSI_JACK, main = 'Jackk', xlim = c(0, 1))
density(GSI_BOOT, main = 'Boot', xlim = c(0, 1))

#par(mfrow = c(1, 1))
#plot(densJack, main = 'GSI values', col = 'blue', ylim = c(0, 10))
#lines(densBoot, col = 'red')
#legend('topleft', legend = c('Jackknife', 'Bootstrap'), 
#       lty = c(1, 1), lwd = c(1, 1), col = c('blue', 2))

# Guardar los resultados en archivos raster
writeRaster(GSI_BOOT,"GSI/INDICE_GSI_BOOT.tif", overwrite = TRUE)
writeRaster(GSI_JACK,"GSI/INDICE_GSI_JACK.tif", overwrite = TRUE)

