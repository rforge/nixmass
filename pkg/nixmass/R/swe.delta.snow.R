swe.delta.snow <- function(data, rho.max=401.2588, rho.null=81.19417, c.ov=0.0005104722, k.ov=0.37856737, k=0.02993175, tau=0.02362476, eta.null=8523356, timestep=24, verbose=FALSE) {
     
  if(!inherits(data,"data.frame"))
    stop("swe.deltasnow: data must be of class data.frame")
  
  if(!any((is.element(colnames(data), c("hs","date")))))
    stop("swe.deltasnow: data must contain at least two columns named 'hs' and 'date'")
  
  Hobs <- data$hs
  if(any(is.na(Hobs)))
    stop("swe.deltasnow: snow depth data must not be NA")
  if(!all(Hobs >= 0))
    stop("swe.deltasnow: snow depth data must not be negative")
  if(!is.numeric(Hobs))
    stop("swe.deltasnow: snow depth data must be numeric")
  if(Hobs[1] != 0)
    stop("swe.deltasnow: snow depth observations must start with 0")
  
  if(!inherits(data$date,"character"))
    stop("swe.deltasnow: date column must be of class character")
  
  z <- zoo(Hobs,as.Date(data$date))
  if(!is.regular(z, strict = TRUE))
    stop("swe.deltasnow: date column must be strictly regular")
  
  snowpack.dd <- NULL
 
  
  H        <- c()                              # modeled total height of snow at any day [m]
  SWE      <- c()                              # modeled total SWE at any day [kg/m2]
  ly       <- 1                                # layer number [-]
  ts       <- timestep * 3600                  # timestep between observations [s]
  g        <- 9.81                             # gravity [ms-2]
  
  # preallocate matrix as days X layers   
  ly.tot   <- length(which(Hobs>0))            # maximum number of layers [-]
  day.tot  <- length(Hobs)                     # total days from first to last snowfall [-]
  h        <- matrix(0,ly.tot,day.tot)         # modeled height of snow in all layers [m]
  swe      <- matrix(0,ly.tot,day.tot)         # modeled swe in all layers [kg/m2]
  age      <- matrix(0,ly.tot,day.tot)         # age of modeled layers [days]

  # helper
  m        <- rep("",day.tot)                  # vector of (verbose) messages
  prec     <- 10^-10                           # precision for arithmetic comparisons [-]
  
 
  
  
  
  #-------------------------------------------------------------------------
  # compaction of snowpack
  #-------------------------------------------------------------------------
  dry_metamorphism <- function(h.d, swe.d, age.d, ly.tot, ly, k, rho.max, ts, prec, g){
    # h.d=h[,t];swe.d=swe[,t];age.d=age[,t]
    
    #snowpack.dd <- NULL
      
    # .d  -> today
    # .dd -> tomorrow
    
    # compute overburden for each layer
    # the overburden for the first layer is the layer itself
    swe.hat.d <- c()
    for(i in 1:ly.tot){
      swe.hat.d[i] <- sum(swe.d[i:ly.tot])
    }
    
    snowpack.d   <- data.frame(h  = h.d, swe = swe.d, swe.hat = swe.hat.d, age = age.d)
    H.d          <- sum(snowpack.d$h)
    
    #snowpack.dd  <<- data.frame(t(apply(snowpack.d, 1, compactH, H.d))) # compress today's snowpack layerwise
    #colnames(snowpack.dd) <<- c("h","swe","age","rho")
    #rownames(snowpack.dd) <<- paste0("dd.layer",1:nrow(snowpack.dd))
    
    a <- data.frame(t(apply(snowpack.d[(1:ly),], 1, compactH, H.d, k, rho.max, ts, prec, g, eta.null)))
    b <- data.frame(rep(0,ly.tot-ly),rep(0,ly.tot-ly),rep(0,ly.tot-ly),rep(0,ly.tot-ly))
    colnames(a) <- colnames(b) <- c("h","swe","age","rho")
    snowpack.dd <<- rbind(a,b)
    rownames(snowpack.dd) <<- paste0("dd.layer",1:nrow(snowpack.dd))
    
    
    return(snowpack.dd)
  }
  
  
  
  #-------------------------------------------------------------------------
  # compaction of individual snow layers without additional load 
  # today's values are compacted to tomorrow's values
  #-------------------------------------------------------------------------
  compactH <- function(x, H.d, k, rho.max, ts, prec, g, eta.null){
    # .d  -> today
    # .dd -> tomorrow
    age.d <- ifelse(x[1] == 0, 0, x[4])
    h.dd <- x[1]/(1 + (x[3] * g * ts)/eta.null * exp(-k * x[2]/x[1]))
    h.dd <- ifelse(x[2]/h.dd > rho.max, x[2]/rho.max, h.dd)
    h.dd <- ifelse(x[1] == 0, 0, h.dd)
    swe.dd  <- x[2]
    age.dd  <- ifelse(x[1] == 0, 0, age.d + 1)
    rho.dd  <- ifelse(x[1] == 0, 0, swe.dd/h.dd)
    rho.dd  <- ifelse(rho.max - rho.dd < prec, rho.max, rho.dd)
    return(cbind(h = h.dd, swe = swe.dd, age = age.dd, rho = rho.dd))
  } 
  
  
  
  #-------------------------------------------------------------------------
  # assigns snowpack properties of tomorrow
  #-------------------------------------------------------------------------
  assignH <- function(sp.dd, h, swe, age, H, SWE, t, day.tot){
    if(t < day.tot){
      h[,t+1]   <- sp.dd$h
      swe[,t+1] <- sp.dd$swe
      age[,t+1] <- sp.dd$age
      
      H[t+1]   <- sum(h[,t+1])
      SWE[t+1] <- sum(swe[,t+1])
    }
    return(list(h=h, swe=swe, age=age, H=H, SWE=SWE))
  }
  
  
  
  drenchH <- function(t, ly, ly.tot, day.tot, Hobs, h, swe, age, H, SWE, rho.max, c.ov, k.ov, k, ts, prec, m){
    Hobs.d=Hobs[t]; h.d=h[,t]; swe.d=swe[,t]; age.d=age[,t]
             
    msg(m,t,paste("melt "))
    
    runoff <- 0
    
    # distribute mass top-down
    for(i in ly:1){
         if( sum(h.d[-i]) + swe.d[i]/rho.max - Hobs.d >= prec ){
              # layers is densified to rho.max
              h.d[i] <- swe.d[i]/rho.max  
         } else {
              # layer is densified as far as possible
              # but doesnt reach rho.max
              h.d[i] <- swe.d[i]/rho.max + abs(sum(h.d[-i]) + swe.d[i]/rho.max - Hobs.d) 
              break
         }
         
    }
    
    # all layers have rho.max
    if( all(rho.max - swe.d[1:ly]/h.d[1:ly] <= prec) ){
      msg(m,t,paste("no further compaction "))
      
      # produce runoff if sum(h.d) - Hobs.d is still > 0
      #if ( sum(h.d) - Hobs.d > prec ){
      msg(m,t,paste("runoff "))
      # decrease swe from all layers?
      # or beginning with lowest?
      #swe.d[1:ly] <- swe.d[1:ly] - (sum(h.d) - Hobs.d) * rho.max
      scale <- Hobs.d/sum(h.d)
      runoff <- (sum(h.d) - Hobs.d) * rho.max  # excess is converted to runoff [kg/m2]
      h.d <- h.d * scale                       # all layers are compressed (and have rho.max) [m]
      swe.d <- swe.d * scale
      
      # }
      
    } else {
      msg(m,t,paste("compaction "))
    }
    
    h[,t]   <- h.d
    swe[,t] <- swe.d
    age[,t] <- age.d
    H[t]    <- sum(h[,t])
    SWE[t]  <- sum(swe[,t])
    
    # no further compaction possible
    #snowpack.tomorrow <- cbind(h = h.d, swe = swe.d, age = age.d, rho = swe.d/h.d)
    #colnames(snowpack.tomorrow) <- c("h","swe","age","rho")
    snowpack.tomorrow <- dry_metamorphism(h[,t], swe[,t], age[,t], ly.tot, ly, k, rho.max, ts, prec, g)
    
    # set values for next day
    rl <- assignH(snowpack.tomorrow, h, swe, age, H, SWE, t, day.tot)
    h <- rl$h
    swe <- rl$swe
    age <- rl$age
    H <- rl$H
    SWE <- rl$SWE
    
    return(list(h=h, swe=swe, age=age, H=H, SWE=SWE))
  }
  
  scaleH <- function(t, ly, ly.tot, day.tot, deltaH, Hobs, h, swe, age, H, SWE, rho.max, k, ts, prec, m){
    
    # re-compact snowpack from yesterdays values with adapted eta
    # .d  -> yesterday
    # .dd -> today
    Hobs.d  = Hobs[t-1]
    Hobs.dd = Hobs[t]
    h.d     = h[,t-1]
    swe.d   = swe[,t-1]
    age.d   = age[,t]  #; deltaH.d = deltaH
    
    # todays overburden   
    swe.hat.d <- c()
    for(i in 1:ly.tot){
      swe.hat.d[i] <- sum(swe.d[i:ly.tot])
    }
    
    # analytical solution for layerwise adapted viskosity eta
    # assumption: recompaction ~ linear height change of yesterdays layers (see paper)
    eta.cor <- c()
    for(i in 1:ly.tot){
         rho.d <- swe.d[i]/h.d[i]
         x <- ts * g * swe.hat.d[i] * exp(-k*rho.d) # yesterday
         P <- h.d[i]/Hobs.d # yesterday
         eta.i <- Hobs.dd * x * P / (h.d[i] - Hobs.dd * P)
         eta.cor <- c(eta.cor, ifelse(is.na(eta.i), 0, eta.i))
    }
    
    # compute H of today with corrected eta
    # so that modeled H = Hobs
    h.dd.cor <- h.d/(1 + (swe.hat.d * g * ts)/eta.cor * exp(-k * swe.d/h.d)) 
    h.dd.cor[which(is.na(h.dd.cor))] <- 0
    H.dd.cor <- sum(h.dd.cor)
    
    # and check, if Hd.cor is the same as Hobs.d
    if(abs(H.dd.cor - Hobs.dd) > prec)
      warning(paste0("day ",t,": error in exponential re-compaction: H.dd.cor-Hobs.dd=",H.dd.cor - Hobs.dd))
   
    
    # which layers exceed rho.max?
    idx.max <- which(swe.d/h.dd.cor - rho.max > prec)
    
    if(length(idx.max) > 0){
         
         if(length(idx.max) < ly){
              # collect excess swe in those layers
              swe.excess <- swe.d[idx.max]-h.dd.cor[idx.max]*rho.max
              
              # set affected layer(s) to rho.max
              swe.d[idx.max] <- swe.d[idx.max] - swe.excess             
              
              # distribute excess swe to other layers top-down
              lys <- 1:ly
              lys <- lys[-idx.max]
              i <- lys[length(lys)]
              swe.excess.all <- sum(swe.excess)
              while(swe.excess.all > 0){
                   swe.res <- h.dd.cor[i] * rho.max - swe.d[i] # layer tolerates this swe amount to reach rho.max
                   if(swe.res > swe.excess.all){
                        swe.res <- swe.excess.all
                   }
                   swe.d[i]  <- swe.d[i] + swe.res
                   swe.excess.all  <- swe.excess.all - swe.res
                   i <- i - 1
                   if(i<=0 & swe.excess.all > 0){
                        msg(m,t,paste(" runoff"))
                        break
                   }
              }
         } else {
              # if all layers have density > rho.max
              # remove swe.excess from all layers (-> runoff)
              # (this sets density to rho.max)
              swe.excess <- swe.d[idx.max]-h.dd.cor[idx.max]*rho.max
              swe.d[idx.max] <- swe.d[idx.max] - swe.excess
              msg(m,t,paste(" runoff"))
         }
    }
    
    # if(any(na.omit(swe.d/h.dd.cor) - rho.max > prec)){
    #      stop()
    # }
    
    h[,t]   <- h.dd.cor
    swe[,t] <- swe.d 
    age[,t] <- age.d
    H[t]    <- sum(h[,t])
    SWE[t]  <- sum(swe[,t])
    
    # compact actual day 
    # if all layers already have maximum density rho.max
    # the snowpack will not be changed by the following step
    snowpack.tomorrow <- dry_metamorphism(h[,t], swe[,t], age[,t], ly.tot, ly, k, rho.max, ts, prec, g)
    
    # set values for next day
    rl <- assignH(snowpack.tomorrow, h, swe, age, H, SWE, t, day.tot)
    h <- rl$h
    swe <- rl$swe
    age <- rl$age
    H <- rl$H
    SWE <- rl$SWE
    
    return(list(h=h, swe=swe, age=age, H=H, SWE=SWE))
    
  }
    
    
  #-------------------------------------------------------------------------
  # keep track of messages in a vector
  #-------------------------------------------------------------------------
  msg <- function(m,t,strg){
    if(verbose){
      cat(paste(strg))
      if(is.null(m[t])){
        m[t] <- strg  
      } else {
        m[t] <- paste(m[t],strg)
      }
    }  
  }
  
  
  
  
  if(verbose){
       cat("Using parameters:\n",
           "rho.max  =",rho.max,"\n",
           "rho.null =",rho.null,"\n",
           "c.ov     =",c.ov,"\n",
           "k.ov     =",k.ov,"\n",
           "k        =",k,"\n",
           "tau      =",tau,"\n",
           "eta.null =",eta.null
           )
  }
  
  for (t in 1:day.tot) {
    
    msg(m,t,paste("day",t,": "))
    
   
    # snowdepth = 0, no snow cover
    if( Hobs[t] == 0 ){      
      if(t > 1){
        if(Hobs[t-1] == 0){
          msg(m, t, paste0(""))        
        } else {
          msg(m, t, paste0("runoff"))
        }        
      }else {
        msg(m, t, paste0(""))         
      }
      H[t]    <- 0
      SWE[t]  <- 0
      h[,t]   <- 0
      swe[,t] <- 0
      
      # there is snow
    } else if( Hobs[t] > 0 ){
      
      # first snow in/during season
      if( Hobs[t-1] == 0 ){
        ly <- 1
        msg(m,t,paste("produce layer",ly))
        age[ly,t]     <- 1
        h[ly,t]       <- Hobs[t]   
        H[t]          <- Hobs[t]
        swe[ly,t]     <- rho.null * Hobs[t]
        SWE[t]        <- swe[ly,t]
        
        
        # compact actual day 
        snowpack.tomorrow <- dry_metamorphism(h[,t], swe[,t], age[,t], ly.tot, ly, k, rho.max, ts, prec, g)
        
        # set values for next day
        rl <- assignH(snowpack.tomorrow, h, swe, age, H, SWE, t, day.tot)
        h <- rl$h
        swe <- rl$swe
        age <- rl$age
        H <- rl$H
        SWE <- rl$SWE
        
        
        # non-first day of snow covered period
      } else if ( Hobs[t-1] > 0 ){
        
        deltaH <- Hobs[t] - H[t]
        
        if( deltaH > tau ){
          msg(m,t,paste("create new layer",ly+1))

          sigma.null <- deltaH * rho.null * g
          epsilon <- c.ov * sigma.null * exp(-k.ov * snowpack.dd$rho/(rho.max - snowpack.dd$rho))
          h[,t]  <- (1 - epsilon) * h[,t]
          #epsilon <- 1 - c.ov * sigma.null * exp(-k.ov * snowpack.dd$rho/(rho.max - snowpack.dd$rho))
          #h[,t]     <- epsilon * h[,t]
          swe[,t]   <- swe[,t-1]
          age[(1:ly),t]   <- age[(1:ly),t-1] + 1
          H[t]      <- sum(h[,t])
          SWE[t]    <- sum(swe[,t])
          #RHO[t]    <- SWE[t]/H[t]
          
          # only for new layer
          ly            <- ly + 1
          h[ly,t]       <- Hobs[t] - H[t]
          swe[ly,t]     <- rho.null * h[ly,t]
          age[ly,t]     <- 1
          
          # recompute
          H[t]   <- sum(h[,t])
          SWE[t] <- sum(swe[,t])
          
          # compact actual day 
          snowpack.tomorrow <- dry_metamorphism(h[,t], swe[,t], age[,t], ly.tot, ly, k, rho.max, ts, prec, g)
          
          # set values for next day
          rl <- assignH(snowpack.tomorrow, h, swe, age, H, SWE, t, day.tot)
          h <- rl$h
          swe <- rl$swe
          age <- rl$age
          H <- rl$H
          SWE <- rl$SWE
          
          
          # no mass gain or loss, but scaling
        } else if( deltaH >= -tau & deltaH <= tau ) {
          msg(m,t,paste("scaling: "))
          rl <- scaleH(t, ly, ly.tot, day.tot, deltaH, Hobs, h, swe, age, H, SWE, rho.max, k, ts, prec, m)
          h <- rl$h
          swe <- rl$swe
          age <- rl$age
          H <- rl$H
          SWE <- rl$SWE
  
          
          
        } else if ( deltaH < -tau ){
          
          msg(m,t,paste("drenching: "))
          rl <- drenchH(t, ly, ly.tot, day.tot, Hobs, h, swe, age, H, SWE, rho.max, c.ov, k.ov, k, ts, prec, m)
          h <- rl$h
          swe <- rl$swe
          age <- rl$age
          H <- rl$H
          SWE <- rl$SWE
          
          
        } else {
          msg(m,t,"??")
        }
        
      }
      
      
    } 
    msg(m,t,"\n")
  }
  
  return(SWE)
  
  
}

  
  
  
  