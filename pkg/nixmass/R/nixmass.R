nixmass <- function(data, model = c("delta.snow","jo09","pi16","st10","gu19"), alt, region.jo09, region.gu19, snowclass.st10, verbose = FALSE) {
   
  model <- match.arg(model, several.ok = TRUE)
  if(length(model) == 0) 
    model <- "deltasnow"

  
  #-----------------------------------------------------------------------
  # split into different models
  swe <- list()
  for(m in model){
   if(m == "delta.snow"){
     swe[["swe"]][[m]] <- swe.delta.snow(data, rho.max=401.2588, rho.null=81.19417, c.ov=0.0005104722, k.ov=0.37856737, k=0.02993175, tau=0.02362476, eta.null=8523356, timestep=24, verbose)
   } else if (m == "jo09"){
     swe[["swe"]][[m]] <- swe.jo09(data, alt, region.jo09)
   } else if (m == "pi16"){
     swe[["swe"]][[m]] <- swe.pi16(data, rho_0=200, K=1)
   } else if (m == "st10"){
     swe[["swe"]][[m]] <- swe.st10(data, snowclass.st10)
   } else if (m == "gu19"){
     swe[["swe"]][[m]] <- swe.gu19(data, region.gu19, n0=NA ,n1=NA, n2=NA)
   } 
  }
  
  swe[["date"]] <- data$date
  swe[["hs"]] <- data$hs
  
  class(swe) <- "nixmass"
  return(swe)
  
}


# S3 function summary
summary.nixmass <- function(object, ...){
  
  if(class(object) != "nixmass")
    stop("nixmass: Object must be of class 'nixmass'.")
  
  models <- names(object$swe)
  if(length(models) == 0)
    stop("nixmass: Cannot plot. No model was computed.")
  
  res <- c()
  for(m in models){
    res <- rbind(res,object$swe[[m]])  
  }
  rownames(res) <- models
  print(apply(res,1,summary))
  
}




# S3 function plot
plot.nixmass <- function(x, title = NULL, ...){
  
  if(class(x) != "nixmass")
    stop("nixmass: Object must be of class 'nixmass'.")
  
  models <- names(x$swe)
  if(length(models) == 0)
    stop("nixmass: Cannot plot. No model was computed.")
  
  colors <- c("#E16A86","#C18500","#799D00","#00AB6E","#00A9BE","#6C8EE6","#D169D0")
  
  Sys.setlocale("LC_TIME", locale = "English")
  
  # define maximum swe for plot outline
  ymax <- c()
  for(m in models){
    ymax <- c(ymax,max(na.omit(c(x$swe[[m]]))))  
  }
  ymax <- max(ymax)
  plot(as.Date(x$date,),xaxt="n", x$swe[[1]],type="n",xlab="",ylab="HS (cm) / SWE (kg/m2)",ylim=c(0,ymax*1.2))
  axis.Date(1, at = seq(as.Date(x$date[1]), as.Date(x$date[length(x$date)]), by = "2 month"), format = "%b")
  n <- 1
  for(m in models){
    lines(as.Date(x$date),x$swe[[m]],type="l",col=colors[n])
    n <- n + 1
  }
  lines(as.Date(x$date),x$hs*100,type="l",lty=2,col="black")
  t <- ifelse(is.null(title),"SWE",title) #paste0("Chartreuse (",alts[s],"m)")
  legend(title=t,"topleft", 
         legend=c(models,"HS"),
         col=c(colors[1:length(models)],"black"), 
         lty=c(rep(1,length(models)),2), cex=0.8, bg="transparent", bty = "n")
  
  invisible(x)
}


