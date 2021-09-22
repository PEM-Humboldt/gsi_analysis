library(raster)
library(fmsb)
library(sf)
library(dismo)
library(rgdal)
library(MASS)
library(ROCR)
library(rgeos)

options(scipen=9999)

#Working Directory
wd <- 'F:/GSI/'
setwd(wd)
dir.create("2_Ambiental/")

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

## remove spp outsite to AOI
system.time(over.coords <- over(Data2, col)); head(over.coords) #531387 records, 8 min
inAOI <- as.numeric(over.coords$OBJECTID)
inAOI[inAOI>=1] <- 1
inAOI[is.na(inAOI)] <- 0
Data2<-as.data.frame(Data2)
Data2$inAOI <- inAOI
Data2<-Data2[Data2$inAOI == 1,]
coordinates(Data2)=~lon+lat
Data2@proj4string@projargs <- GRS.geo

#remove duplicates (originallity to 1 km but to maintain the species records we don´t remove data)
dat_sp1 <- sp::remove.duplicates(Data2, zero = 0.0083333333) # delayed process
crs(dat_sp1)<-GRS.geo
pt_filtered = Data2[0,]       # copy file structure, but no records
spp<-unique(Data2$species)  # store list of species

for(i in 1:length(spp)) {
        print(paste0("Processing ", spp[i], ' (', i, ' to ',length(spp), ' species)'))
        temp  <- subset(Data2, species==spp[i])
        temp1 <-sp::remove.duplicates(temp, zero = 0.0) # to maintan species record dont remove for distance
        pt_filtered = rbind(pt_filtered, temp1)
}
dim(Data)          # 5266539 records (all species records)
dim(pt_filtered) #  366150 records (unique localities for each species)
dim(dat_sp1)   #  231098  records (unique localities)
crs(pt_filtered)<-GRS.geo
pt_filtered2<-as.data.frame(pt_filtered)

###Variables to use
alt1<- raster::getData('worldclim', var="alt", res=.5, lon=-71, lat=2)
alt2<- raster::getData('worldclim', var="alt", res=.5, lon=-71, lat=-4)
alt<-merge(alt1, alt2)
alt<-crop(alt, extent(col))
alt<-mask(alt, col)
plot(alt)
bio1<- raster::getData('worldclim', var="bio", res=.5, lon=-71, lat=2)
bio2<- raster::getData('worldclim', var="bio", res=.5, lon=-71, lat=-4)
bio<-merge(bio1, bio2)
bio<-crop(bio, extent(col))
bio<-mask(bio, col)
rm(alt1, alt2, bio1, bio2)
slope <- terrain(alt, opt='slope', unit='degrees')
var <- stack(bio,alt,slope)
envVarNames <- c(paste0('bio_',seq(1,19)), 'alt', 'slope')
names(var)<-envVarNames
var<-var[[c(1:7,10:17,20,21)]] # Remove bio8,9,18 y 19 for bias (paper XX)

## enviromental aoi with points (to generatin data for model spatial prediction)
buf<- buffer(pt_filtered, width = 500, dissolve=TRUE ) ## buffer 500m
buf<- as(buf,"SpatialPolygonsDataFrame")
gClip <- function(shp, bb){
        if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
        else b_poly <- as(extent(bb), "SpatialPolygons")
        gIntersection(shp, b_poly, byid = TRUE)
}

mas_country<-gClip(buf, col)

#mask to cut but area to analisys
r <- raster(ncol= ncol(var), nrow=nrow(var))
r<-resample(r, var)
mas_rast<- rasterize(mas_country, r)
mas_rast<-resample(mas_rast, colrast, method = 'ngb')
mas_rast <- mask(colrast, mas_rast, inverse = T)

## Random points AOI
random_p<-randomPoints(mas_rast, n= length(dat_sp1))
random <- as.data.frame(random_p)
coordinates(random_p) = ~ x + y
plot(random_p, pch = 16,  cex = 0.00001)
random_p <- SpatialPointsDataFrame(random_p, random)
crs(random_p)<- GRS.geo

### Summatory to pres/pseudoausences 
# Create a 1.5km buffer around the points
pt0_1km = buffer(random_p, width = 1500, dissolve=F)        # buffer on random points
pt1_1km = buffer(pt_filtered, width = 1500, dissolve=F) # buffer on existing points


x<-mask(var, pt0_1km)        
val.pt0 = extract(x=x, y=pt0_1km, fun=mean, df=TRUE)     # extract values for random points
x2<-mask(var, pt1_1km)        
val.pt1 = extract(x=x2, y=pt1_1km, fun=mean, df=TRUE)     # extract values for existing points
val.pt0 = data.frame(na.omit(val.pt0))       # suppress NA values
val.pt1 = data.frame(na.omit(val.pt1))       # idem 
id0 = rep(0, nrow(val.pt0))                  # create id for random points
id1 = rep(1, nrow(val.pt1))                  # create id for existing points
val.pt0 = cbind(id0,val.pt0)                 # combine id and values
val.pt1 = cbind(id1,val.pt1)                 # combine id and values
names(val.pt0)[1] = "sp"                     # change variable name
names(val.pt1)[1] = "sp"                     # idem
val.all = rbind(val.pt0,val.pt1)             # combine data.frames 

group = kfold(val.all, 5)    # create codes for 5 groups
train = val.all[group != 1,] # selecting 80% for training
test  = val.all[group == 1,] # selecting 20% for testing model's accuracity

mod1km = glm(train$sp ~ ., family="binomial", data = train[c(1,3:19)]) # full GLM model
alias(glm(sp ~., family="binomial", data=train[c(1,3:19)]))
var_x=car::vif(mod1km)
mod1km = update(mod1km, . ~ . -bio_6) # dropping variable bio_6 because are aliased in the model
var_x=car::vif(mod1km);var_x
varmod<-vif_func(in_frame = train[c(3:19, 1)], thresh=3, regres= 'lm') #threshold > 3
f_var<-train[,names(train) %in% varmod]
f_var_nms<-colnames(f_var)
cov <-paste((f_var_nms[-1]), collapse = '+')
ec <- (paste(f_var_nms[1], '~', cov))
final_form <- formula(ec)

fmod1km <-glm(final_form, family="binomial", data = f_var) # full GLM model
step<-stepAIC(fmod1km,trace=F, direction = "both")
step$anova
vif(step) ## evsluate if exist other variable to remove >3
summary(step) 

##prediction and performance
new.pred <- predict(fmod1km, test)
pred <-prediction(new.pred,test$sp)
acc.perf<-performance(pred, measure = "acc")
ind<-which.max( slot(acc.perf, "y.values")[[1]] )
# accuracy and threshold to the model
acc<-slot(acc.perf, "y.values")[[1]][ind]
cutoff<-slot(acc.perf, "x.values")[[1]][ind]
print(c(accuracy= acc, cutoff = cutoff))
plot(acc.perf,ylab="Model accuracy", xlab="Cutoff", col="grey50")
abline(v=cutoff, col="red", lty=3)
abline(h=acc, col="red", lty=3)
text(0.70,0.65, paste("Accuracy =", round(acc, 4)), cex=0.7, pos=4)
text(0.77,0.70, paste("Cutoff =", round(cutoff, 4)), cex=0.7, pos=4)
points(cutoff,acc,cex=3)

## Evaluate  presence and ausence
x <- tapply(new.pred, test$sp, mean)
w<-wilcox.test(new.pred ~ test$sp)
boxplot(new.pred ~ test$sp, xlab="Regions", ylab="GLM prediction", col="grey80",
        names=c("w/o records","w/ records"), cex.axis=0.8)
points(x,pch=22,bg="white", cex=1.1)
text(2.1,3.4, paste("W =",w$statistic  , " / p-value < 0.001"), cex=0.8, font=3)

# Model spatialization (for prioritization of survey regions)
modelo  <- predict(var, step, progress='text', index = 1)
plot(modelo)
writeRaster(modelo,"2_Ambiental/ambiental_dimension_GSI.tif", overwrite=TRUE)
rescal_mod<- (modelo - minValue(modelo))/( maxValue(modelo)- minValue(modelo))
rescal_mod<-rescal_mod/maxValue(rescal_mod) # Re-scaling the value using the Max value to obtain values 0 to 1
writeRaster(rescal_mod,"2_Ambiental/ambiental_dimension_GSI_rescal_ajusvar.tif", overwrite=TRUE)
save.image(paste0('2_Ambiental/ambiental_R_object.RData'))
