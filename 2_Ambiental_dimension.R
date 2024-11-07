library(terra)
library(fmsb)
library(sf)
library(dplyr)
library(dismo)
library(MASS)
library(ROCR)
library(vroom)
library(raptr)

# Configuración de opciones
options(scipen=9999)  # Desactiva notación científica para mejorar legibilidad

# Configuración del directorio de trabajo
wd <- '/GSI'
setwd(wd)  # Establece el directorio de trabajo
dir.create("2_Ambiental/")  # Crea un directorio para almacenar los resultados

# Carga de funciones adicionales
source('GAPfunctions.R')

# Definición de sistemas de referencia espacial
GRS.geo <- '+proj=longlat +datum=WGS84 +no_defs'  # Sistema de coordenadas geográficas

# Cargar shapefile del área de estudio
#col <- geodata::gadm(country= "COL", level=0, path= ".")  # función desactivada. Sirve para cargar shapefile de Colombia desde servidor remoto del paquete geodata
col <- read_sf('/Carpeta_shapefile', 'shapefile')
colrast <- rast(col, res=0.008333333)  # Crea un raster con resolución específica
colrast <- terra::rasterize(col, colrast, fun="sum")  # Rasteriza los datos vectoriales


# Carga de registros de especies
Data <- vroom("Archivo_registros.txt", col_names = TRUE)
Data <- Data[, c('gbifID', 'decimalLatitude', 'decimalLongitude', 'species')] #Select ad first column the name of the id record. And for second and third column, latitude and longitude, respectively.
Data <- Data[!is.na(Data$decimalLongitude),]  # Filtra registros sin longitud
Data <- Data[!is.na(Data$decimalLatitude),]   # Filtra registros sin latitud

# Eliminación de duplicados por coordenadas
Data2 <- Data %>% distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Transformación de datos a objeto espacial sf
Data2 <- Data2 %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude')) %>%
  st_set_crs('EPSG:4326')

# Transformación de proyección de coordenadas
Data2 <- sf::st_transform(Data2, crs = GRS.geo)
col=st_as_sf(col)
col <- sf::st_transform(col, crs = GRS.geo) 


# Filtrado de registros dentro del área de interés (AOI)
system.time(over.coords <- st_intersects(Data2, col, sparse = FALSE)); head(over.coords)
Data2 <- Data2[over.coords == TRUE,]

# Eliminación de duplicados dentro de un umbral de distancia
distance_threshold <- 0.0083333333  # Umbral de distancia ajustable
within_distance <- st_is_within_distance(Data2, dist = distance_threshold)
unique_indices <- !duplicated(within_distance)
dat_sp1 <- Data2[unique_indices, ]

# Filtrado de puntos para cada especie
pt_filtered = Data2[0,]       # Copia de la estructura del archivo pero sin registros
spp <- unique(Data2$species)  # Lista de especies únicas

spp <- na.omit(spp)

# Filtrado y eliminación de duplicados por especie
for(i in 1:length(spp)) {
        print(paste0("Processing ", spp[i], ' (', i, ' to ',length(spp), ' species)'))
        temp  <- Data2[Data2$species == spp[i],]
        temp1 <- temp %>% distinct(geometry, .keep_all = TRUE) # to maintan species record dont remove for distance
        pt_filtered = rbind(pt_filtered, temp1)
}

# Resumen de dimensiones de datos
dim(Data)          # Número total de registros
dim(pt_filtered) # Número de localidades únicas por especie
dim(dat_sp1)   # Localidades únicas después de eliminar duplicados


# Transformación de coordenadas de puntos filtrados
pt_filtered <- sf::st_transform(pt_filtered, GRS.geo)
#pt_filtered2<-as.data.frame(pt_filtered)


### Variables ambientales a utilizar
#alt <- geodata::elevation_30s(country= "COL", path= ".") # función desactivada. Sirve para cargar raster de elevación de Colombia desde servidor remoto del paquete geodata
alt=rast('./elevation/COL_elv_msk.tif') # Carga de raster de elevación
alt <- crop(alt, extent(col)) # Recorta la elevación al área de interés
alt <- mask(alt, col) # Aplica máscara de AOI al raster de elevación
plot(alt)

#bio <- geodata::worldclim_country(country= "COL", var='bio', path= ".", version="2.1") # función desactivada. Sirve para cargar capas bioclim (worldclim) para Colombia desde servidor remoto del paquete geodata
bio=rast('./climate/wc2.1_country/COL_wc2.1_30s_bio.tif')
bio<-crop(bio, extent(col))
bio<-mask(bio, col)

slope <- terra::terrain(alt, v='slope', unit='degrees') # Cálculo de pendiente
var <- c(bio,alt,slope) # Combina todas las variables ambientales
envVarNames <- c(paste0('bio_',seq(1,19)), 'alt', 'slope') # Nombres de variables
names(var)<-envVarNames
var<-var[[c(1:7,10:17,20,21)]] # Remueve variables bioclimáticas con sesgo

### Área de interes ambiental (AOI) ambiental con puntos para generar datos para la predicción espacial del modelo

buf<- st_buffer(pt_filtered, dist = 500) ## buffer 500m
buf <- st_as_sf(buf)

#########################
# Intersección del buffer con el área de interés
mas_country <- st_intersection(buf, col)
#########################

# Máscara para recortar el área a analizar
r <- rast(ncol= ncol(var), nrow=nrow(var))  # Crea un raster con las dimensiones del área de estudio
r <- resample(r, var)  # Resmuestrea el raster para que coincida con las variables
mas_rast <- rasterize(mas_country, r, fun="sum")  # Rasteriza la intersección
mas_rast <- resample(mas_rast, colrast, method = 'near')  # Resmuestrea la máscara
mas_rast <- mask(colrast, mas_rast, inverse = TRUE)  # Aplica la máscara inversa

## Generación de puntos aleatorios en el AOI. Estos puntos se usan para el cálculo estadístico de significancia del índice de sesgo.
random_p<-raptr::randomPoints(mas_rast, nrow(dat_sp1)) #Según García Márquez et al. (2012), se deben generar el mismo número de puntos aleatorios al número de localidades de colecta.

random <- as.data.frame(random_p)
random_p <- random %>% st_as_sf(coords = c('x', 'y'), crs = GRS.geo)
plot(random_p, pch = 16,  cex = 0.00001)


### Suma de presencias/pseudoausencias
# Crea un buffer de 1.5 km alrededor de los puntos
pt0_1km <- st_buffer(random_p, dist = 1500)  # Buffer alrededor de los puntos aleatorios
pt1_1km <- st_buffer(pt_filtered, dist = 1500)  # Buffer alrededor de los puntos reales


# Extracción de valores ambientales
val.pt0 <- extract(var, pt0_1km, fun = mean, df = TRUE) # Valores para los puntos aleatorios
val.pt1 <- extract(var, pt1_1km, fun = mean, df = TRUE) # Valores para los puntos reales
val.pt0 <- na.omit(val.pt0)  # Elimina valores NA
val.pt1 <- na.omit(val.pt1)  
id0 <- rep(0, nrow(val.pt0))  # Identificador para los puntos aleatorios
id1 <- rep(1, nrow(val.pt1))  # Identificador para los puntos reales
val.pt0 <- cbind(id0, val.pt0)  # Combina identificador y valores para los puntos aleatorios
val.pt1 <- cbind(id1, val.pt1)  # Combina identificador y valores para los puntos reales
names(val.pt0)[1] <- "sp"  # Renombra columna
names(val.pt1)[1] <- "sp"  
val.all <- rbind(val.pt0, val.pt1)  # Unión de ambos data frames

# División en grupos para validación cruzada
group <- kfold(val.all, k = 5)  # Creación de códigos para cinco categorías
train <- val.all[group != 1,]  # Selecciona el 80% de los datos para entrenamiento
test <- val.all[group == 1,]  # Selecciona el 20% de los datos para probar la precisión del modelo

# Creación del modelo GLM completo
mod1km <- glm(train$sp ~ ., family="binomial", data = train[c(1,3:19)])  # Modelo GLM con todas las variables
alias(glm(sp ~., family="binomial", data=train[c(1,3:19)]))  # Identifica variables colineales

# Evaluación de la multicolinealidad
var_x <- car::vif(mod1km)  # Calcula el Factor de Inflación de la Varianza (VIF)
mod1km <- update(mod1km, . ~ . -bio_6)  # Elimina la variable 'bio_6' por estar aliada en el modelo
var_x <- car::vif(mod1km); var_x  # Recalcula VIF para el modelo actualizado

# Función para seleccionar variables con VIF menor a un umbral
varmod <- vif_func(in_frame = train[c(3:19, 1)], thresh=3, regres= 'lm')  # Umbral de VIF > 3
f_var <- train[,names(train) %in% varmod]  # Filtra las variables seleccionadas
f_var_nms <- colnames(f_var)  # Obtiene los nombres de las variables seleccionadas
cov <- paste((f_var_nms[-1]), collapse = '+')  # Crea una vector con los nombres de las variables
ec <- paste(f_var_nms[1], '~', cov)  # Construcción de la fórmula del modelo
final_form <- formula(ec)  # Convierte el vector en una fórmula

# Modelo GLM final con las variables seleccionadas
fmod1km <- glm(final_form, family="binomial", data = f_var)  # Modelo GLM con las variables seleccionadas
step <- stepAIC(fmod1km, trace=F, direction = "both")  # Selección de variables con el criterio de información de Akaike (AIC)
step$anova  # Muestra el análisis de varianza para los modelos ajustados
vif(step)  # Evalúa si hay otras variables con VIF > 3 para eliminar
summary(step)  # Resumen del modelo final

# Predicción y rendimiento del modelo
new.pred <- predict(fmod1km, test)  # Predice usando el modelo final en el conjunto de prueba
pred <- prediction(new.pred, test$sp)  # Crea un objeto de predicción para evaluar el rendimiento
acc.perf <- performance(pred, measure = "acc")  # Calcula la precisión del modelo
ind <- which.max(slot(acc.perf, "y.values")[[1]])  # Encuentra el índice con la máxima precisión
acc <- slot(acc.perf, "y.values")[[1]][ind]  # Precisión máxima del modelo
cutoff <- slot(acc.perf, "x.values")[[1]][ind]  # Umbral de corte para la predicción
print(c(accuracy= acc, cutoff = cutoff))  # Valores de precisión y el umbral de corte
plot(acc.perf, ylab="Model accuracy", xlab="Cutoff", col="grey50")  # Grafica la precisión del modelo
abline(v=cutoff, col="red", lty=3)  # Punto de umbral de corte
abline(h=acc, col="red", lty=3)  # Punto de precisión máxima
text(0.70,0.65, paste("Accuracy =", round(acc, 4)), cex=0.7, pos=4)  # Valor de precisión
text(0.77,0.70, paste("Cutoff =", round(cutoff, 4)), cex=0.7, pos=4)  # Valor umbral de corte
points(cutoff, acc, cex=3)  # Señala la intersección de precisión máxima y umbral de corte


# Evaluación de la presencia y ausencia
x <- tapply(new.pred, test$sp, mean)  # Cálculo de la media de las predicciones por grupo de presencia/ausencia
w <- wilcox.test(new.pred ~ test$sp)  # Prueba de Wilcoxon para comparar las distribuciones de presencia/ausencia
boxplot(new.pred ~ test$sp, xlab="Regions", ylab="GLM prediction", col="grey80",
        names=c("w/o records", "w/ records"), cex.axis=0.8)  # Boxplot predicciones por grupo
points(x, pch=22, bg="white", cex=1.1)  # Puntos de la media por grupo
text(2.1, 3.4, paste("W =", w$statistic, " / p-value < 0.001"), cex=0.8, font=3)  # Visualización de los estadísticos de Wilcoxon test


# Espacialización del modelo (para priorización de regiones de muestreo)
modelo <- predict(var, step, progress='text', index = 1)  # Predicción espacial del modelo
plot(modelo)  # Visualiza el modelo espacializado
terra::writeRaster(modelo, "2_Ambiental/ambiental_dimension_GSI.tif", overwrite=TRUE)  # Guarda el modelo como un raster

# Reescalado de valores del modelo
rescal_mod <- (modelo - minmax(modelo)[1])/(minmax(modelo)[2] - minmax(modelo)[1])  # Reescalado para valores entre 0 y 1
rescal_mod <- rescal_mod / minmax(rescal_mod)[2]  # Reescalado usando el valor máximo para obtener valores entre 0 y 1
writeRaster(rescal_mod, "2_Ambiental/ambiental_dimension_GSI_rescal_ajusvar.tif", overwrite=TRUE)  # Guarda el modelo reescalado

# Guarda el entorno de trabajo para futuras referencias
save.image(paste0('2_Ambiental/ambiental_R_object.RData'))
