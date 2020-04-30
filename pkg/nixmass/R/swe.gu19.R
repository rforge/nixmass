swe.gu19 <- function(data, region.gu19, n0=NA ,n1=NA, n2=NA){
  
  if(!inherits(data,"data.frame"))
    stop("swe.gu19: data must be given as data.frame")
  
  if(!any((is.element(colnames(data), c("hs","date")))))
    stop("swe.gu19: data must contain at least two columns named 'hs' and 'date'")
  
  Hobs <- data$hs
  if(any(is.na(Hobs)))
    stop("swe.gu19: snow depth data must not be NA")
  if(!all(Hobs >= 0))
    stop("swe.gu19: snow depth data must not be negative")
  if(!is.numeric(Hobs))
    stop("swe.gu19: snow depth data must be numeric")
  
  if(!inherits(data$date,"character"))
    stop("swe.gu19: date column must be of class character")
 
  # check regions 
  if(missing(region.gu19))
       stop("swe.gu19: region.gu19 must be given")
  if (!is.element(region.gu19,c("italy","southwest","central","southeast","myregion")))
       stop("swe.gu19: region.gu19 must be one of 'italy','southwest','central','southeast','myregion'")
  
  guyennet.coeff <- list("italy"     = c(n0=294,   n1=-8.3e-1,  n2=7.7e-3),
                         "southwest" = c(n0=285.9, n1=1.3e-1,   n2=-0.1e-3),
                         "central"   = c(n0=288.9, n1=-9.3e-1,  n2=8.2e-3),
                         "southeast" = c(n0=332.5, n1=-16.8e-1, n2=13.5e-3),
                         "myregion"  = c(n0=n0,    n1=n1,       n2=n2))
  
  # dos: integer day from 1.9. - 31.8.}
  #dos <- ifelse( month(data$date)>8 & month(data$date) <=12, yday(data$date) - 243, yday(data$date) +122 )
  month <- as.numeric( format(as.POSIXct(data$date), "%m"))
  yday  <- as.POSIXlt(data$date, format="%Y-%m-%d")$yday + 1
  dos <- ifelse( month>8 & month<=12, yday-243, yday+122 )
  doy <- dos - 122 # day of year
  
  # doy...days since 1.11.
  d <- data$date
  doys <- c()
  for(i in 1:nrow(data)){
       # if(month(d[i]) %in% c(8,9,10,11,12)){
       #      doy <- as.integer( difftime(d[i],paste0(year(d[i]),"-11-01"), units="days") )
       # } else {
       #      doy <- as.integer( difftime(d[i],paste0(year(d[i])-1,"-11-01"), units="days") )
       # }
       m <- as.numeric( format(as.POSIXct(d[i]), "%m"))
       yd <- as.POSIXlt(d[i], format="%Y-%m-%d")$yday + 1
       if(m %in% c(10,11,12)){
            doy <- yd-366
       } else if (m>=7 & m<=9) { # model isnt able to produce densities in summer 
            doy <- NA
       } else {
            doy <- yd
       }
       
       doys <- c(doys,doy)
  }
  
  if(region.gu19 == "myregion" & any(is.na(c(n0,n1,n2)))){
       stop("swe.gu19: at least one of the coefficients n0, n1, n2 is NULL")
  }
  
  n0  <- guyennet.coeff[[region.gu19]][1]
  n1  <- guyennet.coeff[[region.gu19]][2]
  n2  <- guyennet.coeff[[region.gu19]][3]
  rho <- n0 + n1*(doys + 61) + n2*(doys + 61)^2   # [kg/m3] bulk snow density
  swe <- rho*data$hs                            # [kg/m2] snow water equivalent
  return(swe)
}
