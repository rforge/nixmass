swe.jo09 <- function(data, alt, region.jo09){
  # hs must be in [m]
  
  if(!inherits(data,"data.frame"))
    stop("swe.jo09: data must be given as data.frame")
  
  if(!any((is.element(colnames(data), c("hs","date")))))
    stop("swe.jo09: data must contain at least two columns named 'hs' and 'date'")
  
  Hobs <- data$hs
  if(any(is.na(Hobs)))
    stop("swe.jo09: snow depth data must not be NA")
  if(!all(Hobs >= 0))
    stop("swe.jo09: snow depth data must not be negative")
  if(!is.numeric(Hobs))
    stop("swe.jo09: snow depth data must be numeric")
  
  if(!inherits(data$date,"character"))
    stop("swe.jo09: date column must be of class character")
  
 # z <- zoo(Hobs,as.Date(data$date))
  #if(!is.regular(z, strict = TRUE))
  #  stop("swe.jo09: date column must be strictly regular")
  
  
  #-----------------------------------------------------------------------
  # check alt >= 0
  if(missing(alt))
    stop("swe.jo09: station elevation must be given")
  if (!alt >= 0)
    stop("swe.jo09: station elevation must not be negative")
  if(is.na(alt))
    stop("swe.jo09: station elevation must not be NA")
  if(!is.numeric(alt))
    stop("swe.jo09: station elevation must be numeric")
  
  # check climateclass 1-7
  if(missing(region.jo09))
    stop("swe.jo09: region.jo09 must be given")
  if (!is.element(region.jo09,1:7))
    stop("swe.jo09: region.jo09 must be integer between 1 and 7")
  
  #-----------------------------------------------------------------------
  # Station altitude vs month of observation
  # Coefficients b [kg/m3] and a [kg/m2]. See Jonas et al.(2009), Table 1.
  # "2000" corresponds to altitudes >= 2000m
  # "1400" corresponds to altitudes >= 1400 and < 2000
  # "0"    corresponds to altitudes <  1400
  # "1" to "12" correspond to respective months 
  jonas.coeff<-list(
     "2000"=list("10"=c(NA,NA),"11"=c(206,47),"12"=c(203,52),"01"=c(206,52),"02"=c(217,46),"03"=c(272,26),"04"=c(331,9),"05"=c(378,21),"06"=c(452,8),"07"=c(470,15),"08"=c(NA,NA),"09"=c(NA,NA)),
     "1400"=list("10"=c(NA,NA),"11"=c(183,35),"12"=c(190,47),"01"=c(208,47),"02"=c(218,52),"03"=c(281,31),"04"=c(354,15),"05"=c(409,29),"06"=c(NA,NA),"07"=c(NA,NA),"08"=c(NA,NA),"09"=c(NA,NA)),
     "0"=list("10"=c(NA,NA),"11"=c(149,37),"12"=c(201,26),"01"=c(235,31),"02"=c(279,9),"03"=c(333,3),"04"=c(347,25),"05"=c(413,19),"06"=c(NA,NA),"07"=c(NA,NA),"08"=c(NA,NA),"09"=c(NA,NA))
  )
  
  #----------------------------------------------------------------------- 
  # region.jo09-specific density offset [kg/m3]]. See Jonas et al.(2009), Table 2.
  jonas.offset <- list("0"=0,"1"=7.6,"2"=11.7,"3"=11.8,"4"=-1.1,"5"=-0.3,"6"=12.1,"7"=-14.7)
  #jonas.offset[["8"]] <- mean(unlist(jonas.offset))
    
  jonas.altclass <- function(alt){
    alt <- as.numeric(alt)
    if(alt>=0 & alt<1400){
      return("0")
    } else if (alt>=1400 & alt<2000){
      return("1400")
    } else if (alt>=2000){
      return("2000")
    } else {
      stop(paste("altitude not valid for model: ",alt))
    }
  }

  #-----------------------------------------------------------------------
  calc.swe <- function(x,alt,region.jo09){
    m <- sprintf("%02d",as.numeric( strftime(x[1], format = "%m") ) ) # month
    hs <- as.numeric(x[2])
    c <- jonas.coeff[[jonas.altclass(alt)]][[m]]
    offset <- as.numeric(jonas.offset[as.character(region.jo09)])
    bd <- ifelse( hs==0, NA, c[2] * hs + c[1] + offset) # # kg/m3
    swe <- ifelse( hs==0, 0, bd * hs) # mm or kg/m2
    return(swe)
  }
 
  df <- data.frame(date=data$date,hs=Hobs)
  swe <- apply(df,1,calc.swe,alt,region.jo09)  
  return(swe)
 }
