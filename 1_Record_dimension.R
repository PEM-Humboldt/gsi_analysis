
library(sf)
library(terra)
library(spatstat)
library(geodata)
library(dplyr)
library(vroom)

# ---------------------------
#   Data loading
# ---------------------------
# change
# New chnage

# Definición del directorio de trabajo
wd <- '/Users/elkintenorio/Library/CloudStorage/GoogleDrive-etenorio@humboldt.org.co/Mi\ unidad/Proyectos/Mapa_de_Vacios/Gap_Selection\ Index_GSI'
setwd(wd)  # Establecer el directorio de trabajo
dir.create("1_Records/")  # Crear un directorio para guardar resultados

# Cargar shapefile del área de estudio
##col <- geodata::gadm(country= "COL", level=0, path= ".") # función desactivada. Sirve para cargar shapefile de Colombia desde servidor remoto del paquete geodata
col <- read_sf('/Users/elkintenorio/Library/CloudStorage/GoogleDrive-etenorio@humboldt.org.co/Mi\ unidad/Proyectos/Mapa_de_Vacios/FPV_2024/AMAZONAS', 'AMAZONAS')

# Definición del sistema de referencia geográfico (GCS) y proyectado (CRS)
GRS.geo<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Geographic Reference to AOI
CRS.proj<-'+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' # Planar Reference to AOI

# Transformación del shapefile al sistema de coordenadas proyectadas
proj.col <- sf::st_transform(col, crs = CRS.proj) 

# Convertir el shapefile a una ventana espacial compatible con spatstat
shape_zoneOwin <- as.owin(proj.col)

# Cargar registros de especies desde un archivo txt
Data <- vroom("/Users/elkintenorio/Library/CloudStorage/GoogleDrive-etenorio@humboldt.org.co/Mi\ unidad/Proyectos/Mapa_de_Vacios/FPV_2024/Registros_Agosto_2024/dwc_AMAZONAS.txt", col_names = TRUE)
Data2 <- Data[, c('gbifID', 'decimalLatitude_x', 'decimalLongitude_x')] #Select ad first column the name of the id record. And for second and third column, latitude and longitude, respectively.
Data2 <- Data2[!is.na(Data2$decimalLatitude_x),] # Eliminar registros con latitud NA
Data2 <- Data2[!is.na(Data2$decimalLongitude_x),] # Eliminar registros con longitud NA
# No estoy seguro que lo que aquí se hace sea correcto. Se está eliminando los duplicados por coordenada
# Esto implica que se deja un punto por lcalidad, removiendo todas las especies, y dejando solo una
# Esto puede ser correcto para el segundo código (dimensión ambiental), pero para este (dimensión de registros), no estoy seguro.  


# Eliminar duplicados basados en coordenadas
Data2 <- Data2[!duplicated(Data2[c("decimalLongitude_x", "decimalLatitude_x")]), ]

# Convertir a objeto sf y transformar las coordenadas al sistema proyectado
coords <- st_as_sf(Data2, coords = c("decimalLongitude_x", "decimalLatitude_x"), crs = GRS.geo)
coordinates.col <- st_transform(coords, crs = CRS.proj)

# Intersección espacial para verificar si las coordenadas están dentro del AOI (Area of Interest)
system.time(over.coords <- st_intersects(coordinates.col, proj.col, sparse = FALSE)); head(over.coords) #531387 records, 8 min
Data2 <- Data2[over.coords == TRUE,]
#system.time(over.coords <- st_intersection(coordinates.col, proj.col, sparse = FALSE)); head(over.coords) #Identify coordinates within the AOI
#coordinates.col <- over.coords #Este paso asume que el st_intersect solo genera en su output coordenadas que efextivamente haver overlap (no genera NAs).
# El comando con st_intersection genera el objetvo espacial de los overlap, pero es más demorado.
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
