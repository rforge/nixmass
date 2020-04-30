swe.st10 <- function(data, snowclass.st10 = c("alpine","maritime","prairie","tundra","taiga")){
  
  if(!inherits(data,"data.frame"))
    stop("swe.st10: data must be given as data.frame")
  
  if(!any((is.element(colnames(data), c("hs","date")))))
     stop("swe.st10: data must contain at least two columns named 'hs' and 'date'")
  
  Hobs <- data$hs * 100 # hs must be in [cm]
  if(any(is.na(Hobs)))
    stop("swe.st10: snow depth data must not be NA")
  if(!all(Hobs >= 0))
    stop("swe.st10: snow depth data must not be negative")
  if(!is.numeric(Hobs))
    stop("swe.st10: snow depth data must be numeric")
  
  if(!inherits(data$date,"character"))
    stop("swe.st10: date column must be of class character")
  
  #-----------------------------------------------------------------------
  snowclass.st10 <- match.arg(snowclass.st10)
  if(length(snowclass.st10) == 0)
    stop("swe.st10: snowclass.st10 must be one of 'alpine','maritime','prairie','tundra','taiga'")
  
  sturm.params <- data.frame(
    "den.max" = c(0.5975,0.5979,0.5940,0.3630,0.2170),
    "den0" = c(0.2237,0.2578,0.2332,0.2425,0.2170),
    "k1" = c(0.0012,0.0010,0.0016,0.0029,0.0000),
    "k2" = c(0.0038,0.0038,0.0031,0.0049,0.0000))
  row.names(sturm.params) <- c("alpine","maritime","prairie","tundra","taiga")
  
  d <- as.Date(data$date)
  doys <- c()
  for(i in 1:nrow(data)){
    if(month(d[i]) %in% c(10,11,12)){
      doy <- yday(as.Date(d[i]))-366
    } else if (month(d[i])>=7 & month(d[i])<=9) { # model isnt able to produce densities in summer 
      doy <- NA
    } else {
      doy <- yday(as.Date(d[i]))
    }
    doys <- c(doys,doy)
  }
  
  
  #-----------------------------------------------------------------------
  calc.swe <- function(x,snowclass.st10){
    hs  <- x[1]
    doy <- x[2]
    den.max <- sturm.params[snowclass.st10,]$den.max
    den0    <- sturm.params[snowclass.st10,]$den0
    k1      <- sturm.params[snowclass.st10,]$k1
    k2      <- sturm.params[snowclass.st10,]$k2
    bd      <- ifelse(is.na(doy), NA, (den.max - den0)*(1 - exp(-k1*hs - k2*doy)) + den0) # g/cm3
    swe     <- ifelse( hs==0, 0, bd * hs * 10) # mm or kg/m2                 
    return(swe)
  }
  df <- data.frame(hs=Hobs,doy=doys)
  swe <- apply(df,1,calc.swe,snowclass.st10)
  return(swe) 
}

