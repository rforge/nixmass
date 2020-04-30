swe.pi16 <- function(data, rho_0=200, K=1){
  
  if(!inherits(data,"data.frame"))
    stop("swe.pi16: data must be given as data.frame")
  
  if(!any((is.element(colnames(data), c("hs","date")))))
    stop("swe.pi16: data must contain at least two columns named 'hs' and 'date'")
  
  Hobs <- data$hs
  if(any(is.na(Hobs)))
    stop("swe.pi16: snow depth data must not be NA")
  if(!all(Hobs >= 0))
    stop("swe.pi16: snow depth data must not be negative")
  if(!is.numeric(Hobs))
    stop("swe.pi16: snow depth data must be numeric")
  
  if(!inherits(data$date,"character"))
    stop("swe.pi16: date column must be of class character")
  # 
  # z <- zoo(Hobs,as.Date(data$date))
  # if(!is.regular(z, strict = TRUE))
  #   stop("swe.jonas: date column must be strictly regular")
  # 
  
  
  # dos: integer day from 1.9. - 31.8.}
  #dos <- ifelse( month(data$date)>8 & month(data$date) <=12,yday(data$date) - 243, yday(data$date) +122 )
  month <- as.numeric( format(as.POSIXct(data$date), "%m"))
  yday  <- as.POSIXlt(data$date, format="%Y-%m-%d")$yday + 1
  dos <- ifelse( month>8 & month<=12, yday-243, yday+122 )
  doy <- dos - 122
  # if(month(data$date) %in% c(10,11,12)){
  #      doy <- yday(as.Date(data$date))-366
  # } else if (month(data$date)>=7 & month(data$date)<=9) { # model isnt able to produce densities in summer 
  #      doy <- NA
  # } else {
  #      doy <- yday(data$date)
  # }
  #rho_0 <- 200                 # [kg/m3] bulk snow density at DOY = -62 (31 October)
  #K <- 1                       # [kg/m3/day] rate of (bulk snow) density increase
  rho <- rho_0 + K*(doy + 61)  # [kg/m3] bulk snow density
  swe <- rho*data$hs           # [kg/m2] snow water equivalent
  return(swe)
}
