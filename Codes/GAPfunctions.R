##Function to GSI analysis

#By: Cristian Cruz-Rodr[i]guez / Iv[a]n gonz[a]lez & Elkin Noguera Urbano
#Date: 14-09-2021

#Edited by: Elkin A. Tenorio
#Date: 26-06-2026


# GAP functions

# =========================================================
# Function: vif_func
#
# Used in:
#   - 2_Ambiental_dimension.R
#
# Description:
# This function identifies environmental predictor
# variables with Variance Inflation Factor (VIF) values
# below a selected threshold in order to reduce
# multicollinearity among explanatory variables.
#
# The function iteratively removes variables with the
# highest VIF values until all remaining predictors are
# below the specified threshold.
#
# This procedure is used during GLM calibration for the
# environmental dimension of the GSI.
#
# Inputs:
# - in_frame:
#     Data frame containing predictor variables
#
# - thresh:
#     Maximum acceptable VIF threshold
#
# - regres:
#     Regression type ('lm' or 'glm')
#
# Outputs:
# - A vector containing the names of variables retained
#   after VIF filtering
#
# Main steps:
# 1. Calculate VIF values for all predictors
# 2. Identify the variable with the highest VIF
# 3. Remove highly collinear variables iteratively
# 4. Stop when all variables are below the threshold
# =========================================================

vif_func<-function(in_frame=x, thresh=y, trace=T, regres = r,...){
  require(fmsb)
  require(car)
  
  if(class(in_frame) != 'data.frame') {
    stop ("It's nessesary to include the table as data.frame")
  }
  
  #get initial vif value for all comparisons of variables
  if (regres == 'lm'){
    vif_init<-NULL
    var_names <- names(in_frame)
    for(val in var_names){
      regressors <- var_names[-which(var_names == val)]
      form <- paste((regressors[-length(regressors)]), collapse = '+')
      ec <- (paste(val, '~', form))
      form_in <- formula(ec)
      vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
    }
    vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
    
    if(vif_max < thresh){
      if(trace==T){ #print output of each iteration
        prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
        cat('\n')
        cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
      }
      return(var_names)
    }
    else{
      
      in_dat<-in_frame
      
      #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
      while(vif_max >= thresh){
        
        vif_vals<-NULL
        var_names <- names(in_dat)
        
        for(val in var_names){
          regressors <- var_names[-which(var_names == val)]
          form <- paste((regressors[-length(regressors)]), collapse = '+')
          ec <- (paste(val, '~', form))
          form_in <- formula(ec)
          
          vif_add<-VIF(lm(form_in, data = in_dat))
          vif_vals<-rbind(vif_vals,c(val,vif_add))
        }
        max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
        
        vif_max<-as.numeric(vif_vals[max_row,2])
        
        if(vif_max<thresh) break
        
        if(trace==T){ #print output of each iteration
          prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
          cat('\n')
          cat('ecuation: ', ec, '\n\n')
          cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
          flush.console()
        }
        
        in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
        
      }
      
      return(names(in_dat))
      
    }
    
  } 
  else if  (regres == 'glm'){
    vif_init2<-NULL
    var_names <- names(in_frame)
    for(val in var_names){
      regressors <- var_names[-which(var_names == val)]
      form <- paste((regressors[-length(regressors)]), collapse = '+')
      ec <- (paste(val, '~', form))
      form_in <- formula(ec)
      vif_init2<-rbind(vif_init2, c(val, car::vif(glm(form_in, data = in_frame))))
    }
    vif_max<-max(as.numeric(vif_init2[,2]), na.rm = TRUE)
    
    if(vif_max < thresh){
      if(trace==T){ #print output of each iteration
        prmatrix(vif_init2,collab=c('var','vif'),rowlab=rep('',nrow(vif_init2)),quote=F)
        cat('\n')
        cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
      }
      return(var_names)
    }
    else{
      
      in_dat<-in_frame
      
      #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
      while(vif_max >= thresh){
        
        vif_vals<-NULL
        var_names <- names(in_dat)
        
        for(val in var_names){
          regressors <-var_names[-which(var_names == val)]
          form <-paste((regressors[-length(regressors)]), collapse = '+')
          ec <-(paste(val, '~', form))
          form_in <- formula(ec)
          vif_add <-car::vif(glm(form_in, data = in_dat))
          vif_vals <-rbind(vif_vals,c(val,vif_add))
        }
        max_row<- which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
        
        vif_max<-as.numeric(vif_vals[max_row,2])
        
        if(vif_max<thresh) break
        
        if(trace==T){ #print output of each iteration
          prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
          cat('\n')
          cat(',removed: ',vif_vals[max_row,1],vif_max,'\n\n')
          flush.console()
        }
        
        in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
        
      }
      
      return(names(in_dat))
      
    }
  }
  else {
    stop("you need select 'glm' or 'lm' regression")
  }
}

# =========================================================
# Function: compBoot
#
# Used in:
#   - 3_Complementarity_dimension.R
#
# Description:
# Estimates expected species richness using the
# Bootstrap non-parametric richness estimator.
#
# This function compares observed richness against
# expected richness to evaluate inventory completeness
# within each spatial unit.
#
# Inputs:
# - sppList:
#     Vector containing species records for a sampling unit
#
# Outputs:
# - Estimated species richness using the Bootstrap method
#
# Main steps:
# 1. Calculate observed richness
# 2. Estimate undetected species probability
# 3. Compute Bootstrap richness estimate
# =========================================================

compBoot <- function(sppList){
  Sobs <- length(unique(sppList))
  Sexp <- Sobs + sum((1 - (table(sppList) / length(sppList))) ** length(sppList))
  return(Sexp) 
}


# =========================================================
# Function: compJack
#
# Used in:
#   - 3_Complementarity_dimension.R
#
# Description:
# Estimates species richness using first- or second-order
# Jackknife estimators.
#
# This estimator uses the frequency of rare species
# (especially singletons) to infer the number of
# undetected species in a sampling unit.
#
# Inputs:
# - sppList:
#     Vector of species records
#
# - nSamples:
#     Number of samples or records
#
# - jackOrder:
#     Jackknife order (1 or 2)
#
# Outputs:
# - Estimated species richness using Jackknife estimation
#
# Main steps:
# 1. Calculate observed richness
# 2. Identify singleton species
# 3. Apply Jackknife richness formula
# =========================================================

compJack <- function(sppList, nSamples, jackOrder = 1){
  Sobs <- length(unique(sppList))
  STable <- table(sppList)
  L <- length(which(STable == 1))
  J1 <- L * ((nSamples - 1)/nSamples)
  if (jackOrder == 1){
    Sexp <- Sobs + J1
    return(Sexp) 
  } else  if (jackOrder == 2){
    Sexp <- Sobs + ((L * (nSamples + nSamples - 3))/ nSamples) -
       ((L * (nSamples - 1) ** 2)/(nSamples * (nSamples - 1)))
    return(Sexp)
  }
}

# NOTE:
# This function provides a simple implementation of
# first- and second-order Jackknife richness estimators.
# The current GSI workflow uses SPECIES::jackknife()
# within the richEst() function for richness estimation.

#compJack <- function(sppList, jackOrder){
#  Sobs <- length(unique(spL)) #j1.S
#  sppFreq <- table(spL) #j1.1
#  n <- table(table(data.j))
  #   m <- data.frame(j = names(n), n_j = n[])
  #   n <- apply(m, 2, as.integer)
#  L <- length(sppFreq[sppFreq==1])
#  j1.m=i1
#  jack1=Sobs+L*((j1.m-1)/j1.m)

#}

# =========================================================
# Function: list2Matrix
#
# Used in:
#   - 3_Complementarity_dimension.R
#
# Description:
# Converts a list object into a structured data frame or
# matrix format.
#
# This function is mainly used to organize richness
# estimation outputs into tabular form.
#
# Inputs:
# - inList:
#     Input list object
#
# - nRow / nCol:
#     Desired matrix dimensions
#
# Outputs:
# - Data frame representation of the input list
#
# Main steps:
# 1. Unlist nested elements
# 2. Reshape into matrix format
# 3. Assign row and column names
# =========================================================

list2Matrix <- function(inList, nRow = NULL, nCol = NULL, colNames = NULL, rowNames = NULL){
  if(!is.null(nCol)) {
    outMatrix <- matrix(unlist(inList), ncol = nCol, byrow = TRUE)
  }
  if(!is.null(nRow)) {
    outMatrix <- matrix(unlist(inList), nrow = nRow, byrow = TRUE)
  }
  if(!is.null(colNames)) {colnames(outMatrix) <- colNames}
  if(!is.null(rowNames)) {rownames(outMatrix) <- rowNames}
  outMatrix <- as.data.frame(outMatrix)
}

# =========================================================
# Function: compRar
#
# Used in:
#   - Currently not directly used in the main workflow
#
# Description:
# Performs rarefaction simulations to estimate species
# accumulation curves under standardized sampling effort.
#
# Inputs:
# - sppList:
#     Species records
#
# - simulations:
#     Number of randomizations
#
# - nObs:
#     Number of observations per simulation
#
# Outputs:
# - Mean and variance of rarefied richness accumulation
#
# Main steps:
# 1. Randomly subsample records
# 2. Calculate accumulation curves
# 3. Repeat simulations
# 4. Estimate mean richness accumulation
# =========================================================

#simulations <- 100; nObs <- 50
compRar <- function(sppList, simulations, nObs){
  simMatrix <- matrix(0, ncol = nObs, nrow = simulations)
  for(s in 1:simulations){
    sppRaref <- sppList[sample(1:length(sppList), nObs, replace = FALSE)]
    sppSimRich <- !duplicated(sppRaref)
    acumCurve <- 1:sum(sppSimRich)
    index.j <- which(sppSimRich)
    if(max(acumCurve) != nObs){  
      filledCurve <- fillCurve(acumCurve, index.j, nObs)
      simMatrix[s, ] <- filledCurve
    } else {
      simMatrix[s, ] <- acumCurve
    }
  }
  return(data.frame(Smean = colMeans(simMatrix, na.rm = TRUE), 
                    Svar = apply(simMatrix, 2, function (x) var(x, na.rm = TRUE))
                    )
         )
}


# =========================================================
# Function: fillCurve
#
# Used in:
#   - compRar()
#
# Description:
# Fills missing positions in species accumulation curves
# generated during rarefaction simulations.
#
# Inputs:
# - acum:
#     Species accumulation values
#
# - index:
#     Index positions of observed richness increments
#
# - lenData:
#     Desired curve length
#
# Outputs:
# - Completed accumulation curve
#
# Main steps:
# 1. Identify missing intervals
# 2. Fill gaps using previous accumulation values
# =========================================================
         
#acum <-acumCurve; index <- index.j; lenData <- 25
fillCurve <- function(acum, index, lenData){
  fill <- rep(NA, lenData)
  fill[index] <- acum
  
  compIndex <- data.frame(cbind(index, c(index.j[-1], lenData), acum))
  compIndex <- compIndex[(compIndex[, 1] - compIndex[, 2]) < -1, ]
  if (nrow(compIndex) >0 ){
    for (c in 1:nrow(compIndex)){
      fill[(compIndex[c, 1] + 1):(compIndex[c, 2] - 1)] <- compIndex[c, 3]
    }
  }
  return(fill)
}

# =========================================================
# Function: richEst
#
# Used in:
#   - 3_Complementarity_dimension.R
#
# Description:
# Calculates multiple non-parametric species richness
# estimators for each spatial unit.
#
# The function computes observed richness and several
# richness estimators including:
#
#   - Bootstrap
#   - Jackknife
#   - Chao estimators
#
# This function forms the core of the complementarity
# dimension of the GSI workflow.
#
# Inputs:
# - sppList:
#     Vector containing species identities
#
# - indexID:
#     Vector defining spatial grouping units (cells)
#
# Outputs:
# - Data frame containing richness estimates and
#   associated statistics for each spatial unit
#
# Main steps:
# 1. Group species records by raster cell
# 2. Calculate observed richness
# 3. Estimate richness using multiple estimators
# 4. Organize outputs into tabular format
# =========================================================
         
library(SPECIES)
# Función paraidentificar la riqueza de especies de los registros en la matriz, usando las celdas asignadas
richEst <- function(sppList, indexID){
  sEstimation <- tapply(sppList, INDEX = indexID, 
                  FUN = function(x){
                    #x <- sppList[indexD == indexID[4]]
                    x <- x[ x != '']
                    if(length(x) == 0) x <- 1
                    n <- data.frame(table(table(x)))
                    n[] <- sapply(n, as.integer)
                    colnames(n) <- c('j', 'n_j')

                    J <- tryCatch(jackknife(n, 1), error = function(e)  rep(NA, 5))
                    C1 <- tryCatch(chao1984(n), error = function(e)  rep(NA, 4))
                    C2 <- tryCatch(ChaoLee1992(n), error = function(e)  rep(NA, 8))
                    C3 <- tryCatch(ChaoBunge(n), error = function(e)  rep(NA, 4))
                    
                    return(c(length(x), length(unique(x)), length(x[ x == '']), compBoot(x),
                      unlist(J), unlist(C1), unlist(C2), unlist(C3)))
                  })
  colEst <- c('Nobs', 'Sobs', 'empty', 'Boot', 'JackknifeOrder', 'JNhat', 'JSE', 'JCI1', 'JCI2', 
    'ChaoNhat', 'ChaoSE', 'ChaoCI1', 'Chao1CI2', 'ChaoLNhat1', 'ChaoLNhat2',
    'ChaoLSE1', 'ChaoLSE2', 'ChaoLCI1', 'ChaoLCI2', 'ChaoLCI3', 'ChaoLCI4',
    'ChaoBNhat', 'ChaoBSE', 'ChaoBCI1', 'ChaoBCI2')
  sEstimation <- list2Matrix(sEstimation, nCol = length(colEst),
                       colNames = colEst,
                       rowNames = names(sEstimation)
                       )
  }



# =========================================================
# Function: normalize01
#
# Used in:
#   - Potential auxiliary normalization step for raster
#     layers in the GSI workflow
#
# Description:
# Rescales raster values between 0 and 1 using min-max
# normalization.
#
# Inputs:
# - x:
#     Raster object
#
# - outDir:
#     Optional output path for saving normalized raster
#
# Outputs:
# - Normalized raster with values ranging from 0 to 1
#
# Main steps:
# 1. Extract raster minimum and maximum values
# 2. Apply min-max normalization
# 3. Optionally export raster to disk
# =========================================================
                                   
normalize01 <- function(x, outDir = NULL){
  xNorm <- (x - x@data@min)/(x@data@max - x@data@min)
  if (!is.null(outDir)){
    writeRaster(xNorm, outDir, overwrite=TRUE)
  }  
  return(xNorm)
}
