##Function to GSI analysis

#By: Cristian Cruz-Rodr[i]guez / Iv[a]n gonz[a]lez & Elkin Noguera Urbano
#Date: 14-09-2021

' Necessary functions to obtain the gap analysis
#' 
#' @param x Matriz con las variables y presencias / ausencias
#' @param y Umbral seleccionado para evaluar el VIF en las variables
#' @return Listado con las variables identificadas \code{x} que poseen un VIF inferior al umbral seleccionado
#' @examples

# GAP functions

#Función para identificar las variables que poseen un VIF inferior al umbral seleccionado
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

#Función para calcular los valores de riqueza usa do el estimador de Boopstrap usando la matriz de datos con los registros de las especies
compBoot <- function(sppList){
  Sobs <- length(unique(sppList))
  Sexp <- Sobs + sum((1 - (table(sppList) / length(sppList))) ** length(sppList))
  return(Sexp) 
}
#Función para calcular los valores de riqueza usando el estimador de Jacknife usando la matriz de datos con los registros de las especies
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

compJack <- function(sppList, jackOrder){
  Sobs <- length(unique(spL)) #j1.S
  sppFreq <- table(spL) #j1.1
  n <- table(table(data.j))
  #   m <- data.frame(j = names(n), n_j = n[])
  #   n <- apply(m, 2, as.integer)
  L <- length(sppFreq[sppFreq==1])
  j1.m=i1
  jack1=Sobs+L*((j1.m-1)/j1.m)

}

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

normalize01 <- function(x, outDir = NULL){
  xNorm <- (x - x@data@min)/(x@data@max - x@data@min)
  if (!is.null(outDir)){
    writeRaster(xNorm, outDir, overwrite=TRUE)
  }  
  return(xNorm)
}
